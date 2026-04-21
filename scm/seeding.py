import random
import numpy as np
from abc import ABC, abstractmethod


class SeedingStrategy(ABC):
    """Base class for all seeding strategies.

    Subclasses must implement seed() which returns a list of node indices
    to be initially infected.
    """

    def __init__(self, N, rho_0, links=None, triangles=None):
        self.N = N
        self.rho_0 = rho_0
        self.num_seeds = int(N * rho_0)
        self.links = links
        self.triangles = triangles

    @abstractmethod
    def seed(self) -> list[int]:
        """Return list of initially infected node indices."""
        ...


class RandomSeeding(SeedingStrategy):
    """Seeds initial infection by randomly selecting rho_0 * N nodes."""

    def seed(self) -> list[int]:
        return random.sample(range(self.N), self.num_seeds)


class HighDegreeSeeding(SeedingStrategy):
    """Seeds the nodes with the highest edge degree (most 1-simplex neighbors)."""

    def seed(self) -> list[int]:
        degrees = np.array([len(self.links[i]) for i in range(self.N)])
        return np.argsort(degrees)[-self.num_seeds:][::-1].tolist()


class HighSimplexSeeding(SeedingStrategy):
    """Seeds the nodes with the highest 2-simplex participation count."""

    def seed(self) -> list[int]:
        triangle_counts = np.array([len(self.triangles[i]) for i in range(self.N)])
        return np.argsort(triangle_counts)[-self.num_seeds:][::-1].tolist()


class MultiContagionSeedingStrategy(ABC):
    """Base class for CIC3 seeding strategies.

    seed() returns a list of length C (one seed-set per contagion). All
    returned seed sets must be pairwise disjoint — a node can only be
    initially infected by one contagion.
    """

    def __init__(self, N, num_seeds_per_contagion, links=None, triangles=None):
        self.N = N
        self.num_seeds_per_contagion = list(num_seeds_per_contagion)
        self.C = len(self.num_seeds_per_contagion)
        self.links = links
        self.triangles = triangles

        total = sum(self.num_seeds_per_contagion)
        if total > N:
            raise ValueError(
                f"sum(num_seeds_per_contagion) = {total} exceeds N = {N}; "
                "disjoint seeding requires total seeds <= N"
            )

    @abstractmethod
    def seed(self) -> list[list[int]]:
        """Return one list of seed node indices per contagion."""
        ...

    def _assert_disjoint(self, seed_sets: list[list[int]]) -> None:
        seen = set()
        for c, s in enumerate(seed_sets):
            dup = seen.intersection(s)
            if dup:
                raise ValueError(
                    f"Seed sets not disjoint: contagion {c} shares nodes "
                    f"{sorted(dup)} with an earlier contagion"
                )
            seen.update(s)


class MultiRandomSeeding(MultiContagionSeedingStrategy):
    """Uniformly random disjoint seed sets per contagion."""

    def seed(self) -> list[list[int]]:
        total = sum(self.num_seeds_per_contagion)
        pool = random.sample(range(self.N), total)

        seed_sets = []
        offset = 0
        for k in self.num_seeds_per_contagion:
            seed_sets.append(pool[offset:offset + k])
            offset += k

        self._assert_disjoint(seed_sets)
        return seed_sets


class MultiHighDegreeSeeding(MultiContagionSeedingStrategy):
    """Top edge-degree nodes, round-robin partitioned across contagions.

    Ranks nodes by 1-simplex degree (descending) and hands them out one
    at a time to each contagion in turn until each has its quota filled.
    """

    def seed(self) -> list[list[int]]:
        return _round_robin_partition(
            ranking=_degree_ranking(self.links, self.N),
            num_seeds_per_contagion=self.num_seeds_per_contagion,
            asserter=self._assert_disjoint,
        )


class MultiHighSimplexSeeding(MultiContagionSeedingStrategy):
    """Top 2-simplex participation nodes, round-robin partitioned across
    contagions.
    """

    def seed(self) -> list[list[int]]:
        return _round_robin_partition(
            ranking=_triangle_ranking(self.triangles, self.N),
            num_seeds_per_contagion=self.num_seeds_per_contagion,
            asserter=self._assert_disjoint,
        )


class CommunityLouvainSeeding(MultiContagionSeedingStrategy):
    """Assigns each contagion's seeds to a distinct detected community.

    Pipeline:
      1. Run single-level Louvain on `links` to partition nodes into
         communities (no external labels required).
      2. If Louvain produces fewer than C communities — typical on
         structureless topologies like Erdos-Renyi or BA — silently
         delegate to the same farthest-first + BFS-ball machinery used
         by `CommunityFarthestFirstSeeding`. This keeps the strategy
         well-defined on any graph; on modular graphs Louvain dominates,
         elsewhere it degrades cleanly into the farthest-first baseline.
      3. Otherwise, rank communities descending by internal edge-endpoint
         density (intra_endpoints / total_incident_endpoints). Hub-like
         communities whose members connect mostly outward sink to the
         bottom of the ranking and are dropped when C < num_communities.
      4. Assign contagion c to the rank-c community; within each, seed
         the nodes with highest *internal* degree (edges staying inside
         the detected community) so boundary / hub-adjacent nodes that
         got merged into a small community are deprioritized.
    """

    def __init__(self, N, num_seeds_per_contagion, links=None, triangles=None,
                 louvain_seed=0, max_iter=20):
        super().__init__(N, num_seeds_per_contagion, links, triangles)
        self.louvain_seed = louvain_seed
        self.max_iter = max_iter

    def seed(self) -> list[list[int]]:
        giant = _largest_component(self.links, self.N)
        communities = _louvain_communities(
            self.links, self.N,
            max_iter=self.max_iter, rng_seed=self.louvain_seed,
            candidates=giant,
        )
        if len(communities) < self.C:
            return _farthest_first_seed_sets(
                self.C, self.num_seeds_per_contagion,
                self.links, self.N, self._assert_disjoint,
            )

        ranked = _rank_by_internal_density(communities, self.links)
        internal_deg = _internal_degrees(communities, self.links, self.N)

        # Per-community sorted pool + cursor. Earlier pop-and-discard logic
        # dropped the remainder of any community a contagion spilled out of;
        # cursors preserve that capacity for later contagions. Communities
        # are still consumed in rank order so contagion 0 gets the densest
        # community, contagion 1 the next densest, etc. — the "one contagion
        # per community" intent is preserved when community sizes cooperate;
        # lopsided Louvain partitions (a few huge clusters + many tiny ones)
        # degrade to "fill from densest clusters in rank order" instead of
        # running out of capacity.
        pools = [
            sorted(comm, key=lambda v: -internal_deg[v])
            for comm in ranked
        ]
        pos = [0] * len(pools)
        ci = 0

        seed_sets = []
        for c, k in enumerate(self.num_seeds_per_contagion):
            collected = []
            while len(collected) < k:
                while ci < len(pools) and pos[ci] >= len(pools[ci]):
                    ci += 1
                if ci >= len(pools):
                    raise ValueError(
                        f"Ran out of community capacity seeding contagion "
                        f"{c}: needed {k}, collected {len(collected)}"
                    )
                pool = pools[ci]
                while pos[ci] < len(pool) and len(collected) < k:
                    v = pool[pos[ci]]
                    pos[ci] += 1
                    collected.append(v)
            seed_sets.append(collected)

        self._assert_disjoint(seed_sets)
        return seed_sets


class CommunityFarthestFirstSeeding(MultiContagionSeedingStrategy):
    """Seeds each contagion in a BFS-ball around a maximally-separated center.

    Pipeline:
      1. Pick C centers via farthest-first traversal. The highest-degree
         node serves as a hub-like *anchor* (never a center); the first
         actual center is the node farthest from it, putting contagion 0
         in the periphery rather than at the hub. Subsequent centers each
         maximize the minimum BFS distance to the existing centers.
      2. Around each center, grow a BFS-ball of seeds outward. Within a
         BFS layer, low-degree nodes are preferred over high-degree hub
         nodes, keeping balls tight inside their local community rather
         than drifting through a hub. Forbidden nodes (already claimed by
         earlier contagions) are traversed past but not added, preserving
         disjointness while keeping locality.

    Captures the community-clustering idea directly from topology without
    naming communities. Hub-like nodes are rarely picked as centers because
    their minimum distance to existing centers is low.
    """

    def seed(self) -> list[list[int]]:
        return _farthest_first_seed_sets(
            self.C, self.num_seeds_per_contagion,
            self.links, self.N, self._assert_disjoint,
        )


def _farthest_first_seed_sets(C, num_seeds_per_contagion, links, N, asserter):
    """Shared implementation of farthest-first + BFS-ball seeding.

    Restricts both center selection and BFS-ball expansion to the largest
    connected component of the graph. Nodes outside the giant component
    are typically isolated or in tiny fragments that cannot host a full
    seed set anyway; dropping them also keeps farthest-first from picking
    "infinitely far" fragment nodes as centers (e.g., on Twitter Mutual).
    Used directly by `CommunityFarthestFirstSeeding` and as a fallback by
    `CommunityLouvainSeeding` on graphs where Louvain under-partitions.
    """
    giant = _largest_component(links, N)
    if sum(num_seeds_per_contagion) > len(giant):
        raise ValueError(
            f"Largest connected component has {len(giant)} nodes but "
            f"{sum(num_seeds_per_contagion)} disjoint seeds are required. "
            f"Graph may be too disconnected for this seed budget."
        )
    centers = _farthest_first_centers(C, links, N, candidates=giant)
    giant_set = set(giant)
    seed_sets = []
    claimed = set()
    for c, center in enumerate(centers):
        k = num_seeds_per_contagion[c]
        ball = _bfs_ball(center, k, links, claimed, N)
        if len(ball) < k:
            raise ValueError(
                f"BFS-ball from center {center} for contagion {c} "
                f"yielded only {len(ball)} unclaimed nodes; need {k}"
            )
        seed_sets.append(ball)
        claimed.update(ball)
    asserter(seed_sets)
    return seed_sets


def _largest_component(links, N):
    """Return the list of node IDs in the largest connected component."""
    comp_id = np.full(N, -1, dtype=np.int64)
    largest, largest_size = [], 0
    for start in range(N):
        if comp_id[start] >= 0:
            continue
        cid = start
        comp_id[start] = cid
        comp = [start]
        frontier = [start]
        while frontier:
            next_f = []
            for u in frontier:
                for v in links[u]:
                    if comp_id[v] == -1:
                        comp_id[v] = cid
                        comp.append(v)
                        next_f.append(v)
            frontier = next_f
        if len(comp) > largest_size:
            largest = comp
            largest_size = len(comp)
    return largest


def _degree_ranking(links, N):
    degrees = np.array([len(links[i]) for i in range(N)])
    return np.argsort(degrees)[::-1].tolist()


def _triangle_ranking(triangles, N):
    counts = np.array([len(triangles[i]) for i in range(N)])
    return np.argsort(counts)[::-1].tolist()


def _round_robin_partition(ranking, num_seeds_per_contagion, asserter):
    C = len(num_seeds_per_contagion)
    seed_sets = [[] for _ in range(C)]
    remaining = list(num_seeds_per_contagion)
    i = 0
    for node in ranking:
        if not any(remaining):
            break
        while remaining[i % C] == 0:
            i += 1
        c = i % C
        seed_sets[c].append(node)
        remaining[c] -= 1
        i += 1
    asserter(seed_sets)
    return seed_sets


def _louvain_communities(links, N, max_iter=20, rng_seed=0, candidates=None):
    """Single-level Louvain modularity optimization (unweighted).

    Δ-modularity for adding isolated node i to community c:
        Δ ∝ k_{i,c} − k_i · Σ_tot[c] / (2m)
    where k_{i,c} = edges from i into c, Σ_tot[c] = sum of degrees of
    nodes currently in c (with i removed), m = |E|.

    If `candidates` is given (typically the largest connected component),
    only those nodes participate in modularity moves and show up in the
    returned communities; non-candidates are filtered out at the end so
    callers don't have to account for singleton fragments.

    Returns: list of communities, each a list of node indices.
    """
    rng = np.random.default_rng(rng_seed)
    if candidates is None:
        cand_mask = np.ones(N, dtype=bool)
    else:
        cand_mask = np.zeros(N, dtype=bool)
        for v in candidates:
            cand_mask[v] = True

    degrees = np.array([len(links[i]) for i in range(N)], dtype=float)
    two_m = degrees[cand_mask].sum()
    if two_m == 0:
        return [[i] for i in range(N) if cand_mask[i]]

    labels = np.arange(N)
    sigma_tot = np.where(cand_mask, degrees, 0.0)
    cand_idx = np.where(cand_mask)[0]

    for _ in range(max_iter):
        changed = False
        for idx in rng.permutation(cand_idx):
            i = int(idx)
            if degrees[i] == 0:
                continue
            c_old = int(labels[i])
            k_i = degrees[i]

            k_in = {}
            for j in links[i]:
                if not cand_mask[j]:
                    continue
                cj = int(labels[j])
                k_in[cj] = k_in.get(cj, 0) + 1

            sigma_tot[c_old] -= k_i
            if c_old not in k_in:
                k_in[c_old] = 0

            best_c, best_gain = c_old, -np.inf
            for c, k_ic in k_in.items():
                gain = k_ic - k_i * sigma_tot[c] / two_m
                if gain > best_gain:
                    best_gain = gain
                    best_c = c

            sigma_tot[best_c] += k_i
            if best_c != c_old:
                labels[i] = best_c
                changed = True
        if not changed:
            break

    groups = {}
    for i in cand_idx:
        groups.setdefault(int(labels[i]), []).append(int(i))
    return list(groups.values())


def _internal_degrees(communities, links, N):
    """Per-node count of neighbors in the same detected community.

    Inside-community connectivity is a better "anchor for this cluster"
    score than raw degree — hub-adjacent boundary nodes that got merged
    into a small community have low internal degree despite high total
    degree, so ranking by internal degree prefers true community cores.
    """
    node_to_comm = np.full(N, -1, dtype=np.int64)
    for idx, comm in enumerate(communities):
        for v in comm:
            node_to_comm[v] = idx
    internal = np.zeros(N, dtype=np.int64)
    for v in range(N):
        cv = node_to_comm[v]
        if cv < 0:
            continue
        for u in links[v]:
            if node_to_comm[u] == cv:
                internal[v] += 1
    return internal


def _rank_by_internal_density(communities, links):
    """Sort communities descending by intra_endpoints / total_endpoints.

    Each undirected intra-community edge contributes 2 to intra_endpoints
    (once per endpoint), matching the sum-of-degrees denominator. Hub-like
    communities whose members mostly connect outward score low.
    """
    node_to_comm = {}
    for idx, comm in enumerate(communities):
        for v in comm:
            node_to_comm[v] = idx

    scored = []
    for comm in communities:
        cid = node_to_comm[comm[0]]
        intra = 0
        total = 0
        for v in comm:
            total += len(links[v])
            for u in links[v]:
                if node_to_comm.get(u) == cid:
                    intra += 1
        density = intra / total if total > 0 else 0.0
        scored.append((density, comm))
    scored.sort(key=lambda item: -item[0])
    return [c for _, c in scored]


def _bfs_distances(start, links, N):
    """BFS from `start`; return per-node distance array (unreachable = -1)."""
    dist = np.full(N, -1, dtype=np.int64)
    dist[start] = 0
    frontier = [start]
    while frontier:
        next_frontier = []
        for u in frontier:
            for v in links[u]:
                if dist[v] == -1:
                    dist[v] = dist[u] + 1
                    next_frontier.append(v)
        frontier = next_frontier
    return dist


def _farthest_first_centers(C, links, N, candidates=None):
    """Pick up to C centers via farthest-first traversal.

    `candidates` (iterable of node IDs or None): restricts center selection
    to these nodes. Callers pass a single connected component so every
    candidate is reachable from every other, keeping BFS distances finite.
    When None (default), all N nodes are eligible — safe only for connected
    graphs.

    Anchor (not a center): the highest-degree candidate — a proxy for the
    cluster's connectivity core. In hub-heavy SBMs this is typically a hub
    node, which we specifically want to avoid using as a contagion center.
    First actual center: argmax BFS distance from the anchor, so contagion
    0 starts in the periphery rather than at the hub.
    Each subsequent center: argmax over unpicked candidates of min BFS
    distance to existing centers; degree breaks ties. Non-candidates and
    unreachable nodes are sentineled at -1 and never picked — unlike a
    naive +inf treatment that would select nodes in disconnected fragments.
    """
    degrees = np.array([len(links[i]) for i in range(N)])
    if candidates is None:
        cand_mask = np.ones(N, dtype=bool)
    else:
        cand_mask = np.zeros(N, dtype=bool)
        for v in candidates:
            cand_mask[v] = True

    cand_idx = np.where(cand_mask)[0]
    if len(cand_idx) == 0 or C <= 0:
        return []

    anchor = int(cand_idx[np.argmax(degrees[cand_idx])])
    d_anchor = _bfs_distances(anchor, links, N)
    # Candidates reachable from the anchor and distinct from it.
    valid = cand_mask & (d_anchor > 0)
    if not valid.any():
        return [anchor][:C]

    max_d = d_anchor[valid].max()
    tied = np.where(valid & (d_anchor == max_d))[0]
    first = int(tied[np.argmax(degrees[tied])]) if len(tied) > 1 else int(tied[0])
    centers = [first]

    d0 = _bfs_distances(first, links, N).astype(float)
    NEG = -1.0  # sentinel: non-candidate, unreachable, or already picked
    min_dist = np.where(cand_mask & (d0 >= 0), d0, NEG)
    min_dist[first] = NEG

    for _ in range(C - 1):
        if min_dist.max() <= 0:
            break
        max_d = min_dist.max()
        tied = np.where(min_dist == max_d)[0]
        nxt = int(tied[np.argmax(degrees[tied])]) if len(tied) > 1 else int(tied[0])
        centers.append(nxt)
        new_d = _bfs_distances(nxt, links, N).astype(float)
        mask = (min_dist > NEG) & (new_d >= 0)
        min_dist[mask] = np.minimum(min_dist[mask], new_d[mask])
        min_dist[nxt] = NEG

    return centers


def _bfs_ball(center, num_seeds, links, forbidden, N):
    """BFS outward from `center`; collect first `num_seeds` unforbidden nodes.

    Within each BFS layer, nodes are considered in ascending-degree order
    so low-degree community members are preferred over high-degree hub
    nodes that happen to sit at the same distance. This keeps balls tight
    within a community rather than drifting through a hub on the first
    layer.

    Forbidden nodes (already claimed by earlier contagions) are traversed
    but not added to the ball, so disjointness is maintained while
    locality around the center is preserved.
    """
    degrees = [len(links[i]) for i in range(N)]
    ball = []
    visited = np.zeros(N, dtype=bool)
    visited[center] = True
    if center not in forbidden:
        ball.append(center)
        if len(ball) == num_seeds:
            return ball
    frontier = [center]
    while frontier:
        layer = []
        for u in frontier:
            for v in links[u]:
                if visited[v]:
                    continue
                visited[v] = True
                layer.append(v)
        layer.sort(key=lambda v: (degrees[v], v))
        for v in layer:
            if v not in forbidden:
                ball.append(v)
                if len(ball) == num_seeds:
                    return ball
        frontier = layer
    return ball
