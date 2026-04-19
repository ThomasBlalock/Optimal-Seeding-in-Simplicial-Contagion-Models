import gzip
import math
import os
import random
import urllib.request

import numpy as np


def _ensure_edge(links, added_edges, i, j):
    """Add edge (i,j) to links/added_edges if not present. Idempotent."""
    edge = (i, j) if i < j else (j, i)
    if edge not in added_edges:
        added_edges.add(edge)
        links[i].append(j)
        links[j].append(i)
        return True
    return False


def _add_triangle(triangles, added_triangles, links, added_edges, i, j, k):
    """Register a 2-simplex on nodes (i,j,k) and ensure its edges exist.

    Returns True if the triangle was newly added.
    """
    nodes = tuple(sorted((i, j, k)))
    if nodes in added_triangles:
        return False
    added_triangles.add(nodes)
    a, b, c = nodes
    triangles[a].append((b, c))
    triangles[b].append((a, c))
    triangles[c].append((a, b))
    for u, v in ((a, b), (b, c), (a, c)):
        _ensure_edge(links, added_edges, u, v)
    return True


def _realized_degrees(N, links, triangles):
    k_avg = sum(len(l) for l in links) / N
    k_delta_avg = sum(len(t) for t in triangles) / N
    return k_avg, k_delta_avg


class RSCGenerator:
    """Generates random simplicial complexes (Erdos-Renyi-style).

    Takes k_avg (mean degree from edges) and k_delta_avg (mean degree from
    triangles). Computes p_1 and p_delta probabilities, then samples edges and
    triangles from binomial distributions. Triangles ensure their constituent
    edges exist (added to the link list if missing).
    """

    def __init__(self, k_avg, k_delta_avg, N=2000):
        self.N = N
        self.k_avg_target = k_avg
        self.k_delta_avg_target = k_delta_avg

        self.p_1 = (k_avg - 2 * k_delta_avg) / (N - 1 - 2 * k_delta_avg)
        self.p_delta = (2 * k_delta_avg) / ((N - 1) * (N - 2))

        self.links = [[] for _ in range(N)]
        self.triangles = [[] for _ in range(N)]
        self.added_edges = set()
        self.added_triangles = set()

    def generate(self, seed=None):
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        self._generate_1_simplices()
        self._generate_2_simplices()

        self.k_avg, self.k_delta_avg = _realized_degrees(
            self.N, self.links, self.triangles
        )
        print(
            f"Realized k_avg = {self.k_avg:.2f}, "
            f"k_delta_avg = {self.k_delta_avg:.2f}"
        )
        return self.links, self.triangles

    def _generate_1_simplices(self):
        print(f"Sampling edges with p_1 = {self.p_1:.8f}")
        num_edges = np.random.binomial(math.comb(self.N, 2), self.p_1)

        node_list = range(self.N)
        while len(self.added_edges) < num_edges:
            i, j = random.sample(node_list, 2)
            _ensure_edge(self.links, self.added_edges, i, j)
            print(
                f"\rEdges sampled: {len(self.added_edges)}/{num_edges}", end=""
            )
        print()

    def _generate_2_simplices(self):
        print(f"Sampling triangles with p_delta = {self.p_delta:.8f}")
        num_triangles = np.random.binomial(math.comb(self.N, 3), self.p_delta)

        node_list = range(self.N)
        while len(self.added_triangles) < num_triangles:
            i, j, k = random.sample(node_list, 3)
            _add_triangle(
                self.triangles, self.added_triangles,
                self.links, self.added_edges, i, j, k,
            )
            print(
                f"\rTriangles sampled: "
                f"{len(self.added_triangles)}/{num_triangles}",
                end="",
            )
        print()


class BAGenerator:
    """Barabasi-Albert preferential attachment with simplicial extension.

    Each arriving node v connects to `m` existing nodes with probability
    proportional to current 1-simplex degree (standard BA). Additionally, v
    initiates Poisson(m_delta) 2-simplices; for each, two existing nodes are
    sampled jointly with PA-biased probability (without replacement) to form
    a triangle (v, i, j). Any missing constituent edges are added to the
    link list (RSC convention), so k_delta_avg inflates k_avg slightly.

    The PA bias uses total 1-simplex degree (including edges added via
    triangles), so hub structure compounds across both simplex orders.
    """

    def __init__(self, m, m_delta, N=2000):
        if m < 1:
            raise ValueError("m must be >= 1")
        if m_delta < 0:
            raise ValueError("m_delta must be >= 0")
        self.N = N
        self.m = m
        self.m_delta = m_delta

        self.links = [[] for _ in range(N)]
        self.triangles = [[] for _ in range(N)]
        self.added_edges = set()
        self.added_triangles = set()

    def generate(self, seed=None):
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        seed_size = self.m + 1
        for i in range(seed_size):
            for j in range(i + 1, seed_size):
                _ensure_edge(self.links, self.added_edges, i, j)

        for v in range(seed_size, self.N):
            self._attach_edges(v)
            self._attach_triangles(v)
            print(f"\rBA nodes added: {v + 1}/{self.N}", end="")
        print()

        self.k_avg, self.k_delta_avg = _realized_degrees(
            self.N, self.links, self.triangles
        )
        print(
            f"Realized k_avg = {self.k_avg:.2f}, "
            f"k_delta_avg = {self.k_delta_avg:.2f}"
        )
        return self.links, self.triangles

    def _pa_sample(self, v, k, exclude=()):
        """Sample k distinct existing nodes (index < v) with probability
        proportional to current 1-simplex degree. Returns [] if infeasible.
        """
        available = v - len(exclude)
        if available < k:
            return []
        degrees = np.fromiter(
            (len(self.links[u]) for u in range(v)),
            dtype=float, count=v,
        )
        for u in exclude:
            if u < v:
                degrees[u] = 0.0
        total = degrees.sum()
        if total <= 0:
            candidates = [u for u in range(v) if u not in exclude]
            if len(candidates) < k:
                return []
            return random.sample(candidates, k)
        probs = degrees / total
        idx = np.random.choice(v, size=k, replace=False, p=probs)
        return idx.tolist()

    def _attach_edges(self, v):
        for u in self._pa_sample(v, self.m):
            _ensure_edge(self.links, self.added_edges, v, u)

    def _attach_triangles(self, v):
        n_tri = np.random.poisson(self.m_delta)
        for _ in range(n_tri):
            partners = self._pa_sample(v, 2)
            if len(partners) < 2:
                continue
            i, j = partners
            _add_triangle(
                self.triangles, self.added_triangles,
                self.links, self.added_edges, v, i, j,
            )


class SBMGenerator:
    """Stochastic Block Model with simplicial extension.

    Edges are sampled per block-pair (a,b) at probability `block_matrix[a,b]`.
    Triangles are sampled over unique block-triples (a<=b<=c) at probability
    `triangle_block_tensor[a,b,c]`. As a shorthand, `triangle_block_probs`
    (length K) specifies intra-community triangle probabilities only.

    Triangles add any missing constituent edges (RSC convention).
    """

    def __init__(self, community_sizes, block_matrix,
                 triangle_block_tensor=None, triangle_block_probs=None):
        self.community_sizes = list(community_sizes)
        self.N = sum(self.community_sizes)
        self.K = len(self.community_sizes)

        self.block_matrix = np.asarray(block_matrix, dtype=float)
        if self.block_matrix.shape != (self.K, self.K):
            raise ValueError(
                f"block_matrix must have shape ({self.K},{self.K})"
            )

        if (triangle_block_tensor is not None) == (
            triangle_block_probs is not None
        ):
            raise ValueError(
                "Provide exactly one of triangle_block_tensor or "
                "triangle_block_probs"
            )

        if triangle_block_tensor is not None:
            t = np.asarray(triangle_block_tensor, dtype=float)
            if t.shape != (self.K, self.K, self.K):
                raise ValueError(
                    f"triangle_block_tensor must have shape "
                    f"({self.K},{self.K},{self.K})"
                )
            self.triangle_block_tensor = t
        else:
            probs = np.asarray(triangle_block_probs, dtype=float)
            if probs.shape != (self.K,):
                raise ValueError(
                    f"triangle_block_probs must have length {self.K}"
                )
            t = np.zeros((self.K, self.K, self.K))
            for a in range(self.K):
                t[a, a, a] = probs[a]
            self.triangle_block_tensor = t

        self.community_nodes = []
        offset = 0
        for size in self.community_sizes:
            self.community_nodes.append(list(range(offset, offset + size)))
            offset += size

        self.links = [[] for _ in range(self.N)]
        self.triangles = [[] for _ in range(self.N)]
        self.added_edges = set()
        self.added_triangles = set()

    def generate(self, seed=None):
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        self._generate_1_simplices()
        self._generate_2_simplices()

        self.k_avg, self.k_delta_avg = _realized_degrees(
            self.N, self.links, self.triangles
        )
        print(
            f"Realized k_avg = {self.k_avg:.2f}, "
            f"k_delta_avg = {self.k_delta_avg:.2f}"
        )
        return self.links, self.triangles

    def _generate_1_simplices(self):
        print("Sampling edges per block-pair")
        total = 0
        for a in range(self.K):
            for b in range(a, self.K):
                p = self.block_matrix[a, b]
                if p <= 0:
                    continue
                nodes_a = self.community_nodes[a]
                nodes_b = self.community_nodes[b]
                if a == b:
                    n_pairs = math.comb(len(nodes_a), 2)
                else:
                    n_pairs = len(nodes_a) * len(nodes_b)
                n_edges = int(np.random.binomial(n_pairs, p))
                added = 0
                while added < n_edges:
                    if a == b:
                        i, j = random.sample(nodes_a, 2)
                    else:
                        i = random.choice(nodes_a)
                        j = random.choice(nodes_b)
                    if _ensure_edge(self.links, self.added_edges, i, j):
                        added += 1
                total += added
                print(
                    f"\rSBM edges sampled (through block {a},{b}): {total}",
                    end="",
                )
        print()

    def _generate_2_simplices(self):
        print("Sampling triangles per block-triple")
        total = 0
        for a in range(self.K):
            for b in range(a, self.K):
                for c in range(b, self.K):
                    p = self.triangle_block_tensor[a, b, c]
                    if p <= 0:
                        continue
                    n_triples = self._count_triples(a, b, c)
                    n_tri = int(np.random.binomial(n_triples, p))
                    added = 0
                    attempts = 0
                    max_attempts = max(10 * n_tri, 1000)
                    while added < n_tri and attempts < max_attempts:
                        triple = self._sample_triple(a, b, c)
                        attempts += 1
                        if triple is None:
                            break
                        i, j, k = triple
                        if _add_triangle(
                            self.triangles, self.added_triangles,
                            self.links, self.added_edges, i, j, k,
                        ):
                            added += 1
                    total += added
                    print(
                        f"\rSBM triangles sampled "
                        f"(through block {a},{b},{c}): {total}",
                        end="",
                    )
        print()

    def _count_triples(self, a, b, c):
        sa = len(self.community_nodes[a])
        sb = len(self.community_nodes[b])
        sc = len(self.community_nodes[c])
        if a == b == c:
            return math.comb(sa, 3)
        if a == b:
            return math.comb(sa, 2) * sc
        if b == c:
            return sa * math.comb(sb, 2)
        return sa * sb * sc

    def _sample_triple(self, a, b, c):
        na = self.community_nodes[a]
        nb = self.community_nodes[b]
        nc = self.community_nodes[c]
        if a == b == c:
            if len(na) < 3:
                return None
            return tuple(random.sample(na, 3))
        if a == b:
            if len(na) < 2 or len(nc) < 1:
                return None
            i, j = random.sample(na, 2)
            return (i, j, random.choice(nc))
        if b == c:
            if len(na) < 1 or len(nb) < 2:
                return None
            j, k = random.sample(nb, 2)
            return (random.choice(na), j, k)
        if not (na and nb and nc):
            return None
        return (random.choice(na), random.choice(nb), random.choice(nc))


class TwitterEgoGenerator:
    """Loads the SNAP ego-Twitter combined edge list and promotes a random
    subset of its triangles to 2-simplices.

    Downloads ``twitter_combined.txt.gz`` from SNAP on first use and caches
    the raw archive under ``data_dir``. Node IDs are relabeled to 0..N-1
    contiguous integers; N is determined by the data (not a constructor
    arg).

    ``edge_mode`` controls directed -> undirected conversion:
    - ``"collapse"``: any directed edge u->v or v->u yields an undirected
      edge {u, v}.
    - ``"mutual"``: {u, v} is kept only if both u->v and v->u are present
      in the SNAP data.

    After ``links`` is built, every triangle in the resulting undirected
    graph is enumerated and independently promoted to a 2-simplex with
    probability ``p_promote``. Non-promoted triangles remain present as
    the three underlying edges only (no entry in ``triangles``).
    """

    DEFAULT_URL = "https://snap.stanford.edu/data/twitter_combined.txt.gz"
    DEFAULT_FILENAME = "twitter_combined.txt.gz"

    def __init__(self, p_promote, edge_mode="collapse",
                 data_dir="data", url=None, filename=None):
        if not 0.0 <= p_promote <= 1.0:
            raise ValueError("p_promote must be in [0, 1]")
        if edge_mode not in ("collapse", "mutual"):
            raise ValueError("edge_mode must be 'collapse' or 'mutual'")
        self.p_promote = p_promote
        self.edge_mode = edge_mode
        self.data_dir = data_dir
        self.url = url or self.DEFAULT_URL
        self.filename = filename or self.DEFAULT_FILENAME

        self.N = None
        self.links = None
        self.triangles = None
        self.added_edges = None
        self.added_triangles = None
        self.k_avg = None
        self.k_delta_avg = None

    def generate(self, seed=None):
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        path = self._fetch_cached()
        edges = self._parse_edges(path)
        self._build_links(edges)
        self._enumerate_and_promote()

        self.k_avg, self.k_delta_avg = _realized_degrees(
            self.N, self.links, self.triangles
        )
        print(
            f"Realized k_avg = {self.k_avg:.2f}, "
            f"k_delta_avg = {self.k_delta_avg:.2f}"
        )
        return self.links, self.triangles

    def _fetch_cached(self):
        os.makedirs(self.data_dir, exist_ok=True)
        path = os.path.join(self.data_dir, self.filename)
        if os.path.exists(path):
            print(f"Using cached Twitter edge list at {path}")
        else:
            print(f"Downloading {self.url} -> {path}")
            urllib.request.urlretrieve(self.url, path)
        return path

    def _parse_edges(self, path):
        print(f"Parsing {path} (edge_mode={self.edge_mode})")
        node_ids = {}

        def intern(raw):
            idx = node_ids.get(raw)
            if idx is None:
                idx = len(node_ids)
                node_ids[raw] = idx
            return idx

        if self.edge_mode == "collapse":
            edges = set()
            with gzip.open(path, "rt") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) < 2:
                        continue
                    u = intern(parts[0])
                    v = intern(parts[1])
                    if u == v:
                        continue
                    edges.add((u, v) if u < v else (v, u))
        else:
            seen = {}
            with gzip.open(path, "rt") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) < 2:
                        continue
                    u = intern(parts[0])
                    v = intern(parts[1])
                    if u == v:
                        continue
                    if u < v:
                        key, bit = (u, v), 1
                    else:
                        key, bit = (v, u), 2
                    seen[key] = seen.get(key, 0) | bit
            edges = {e for e, mask in seen.items() if mask == 3}

        self.N = len(node_ids)
        print(
            f"Parsed {self.N} nodes, {len(edges)} undirected edges "
            f"(edge_mode={self.edge_mode})"
        )
        return edges

    def _build_links(self, edges):
        self.links = [[] for _ in range(self.N)]
        self.triangles = [[] for _ in range(self.N)]
        self.added_edges = set()
        self.added_triangles = set()
        for a, b in edges:
            _ensure_edge(self.links, self.added_edges, a, b)

    def _enumerate_and_promote(self):
        print(f"Enumerating triangles (p_promote={self.p_promote})")
        neighbors = [set(lst) for lst in self.links]
        edges_snapshot = list(self.added_edges)
        found = 0
        promoted = 0
        for u, v in edges_snapshot:
            nu, nv = neighbors[u], neighbors[v]
            if len(nu) <= len(nv):
                small, large = nu, nv
            else:
                small, large = nv, nu
            for w in small:
                if w > v and w in large:
                    found += 1
                    if random.random() < self.p_promote:
                        if _add_triangle(
                            self.triangles, self.added_triangles,
                            self.links, self.added_edges, u, v, w,
                        ):
                            promoted += 1
                    if found % 200000 == 0:
                        print(
                            f"\rTriangles found: {found}  "
                            f"promoted: {promoted}",
                            end="",
                        )
        print(
            f"\rTriangles found: {found}  promoted: {promoted}"
            + " " * 20
        )
