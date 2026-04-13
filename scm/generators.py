import numpy as np
import random
import math


class RSCGenerator:
    """Generates random simplicial complexes (Erdos-Renyi-style).

    Takes k_avg (mean degree from edges) and k_delta_avg (mean degree from
    triangles). Computes p_1 and p_delta probabilities, then samples edges and
    triangles from binomial distributions.
    """

    def __init__(self, k_avg, k_delta_avg, N=2000):
        self.N = N
        self.k_avg = k_avg
        self.k_delta_avg = k_delta_avg

        self.p_1 = (k_avg - 2 * k_delta_avg) / (N - 1 - 2 * k_delta_avg)
        self.p_delta = (2 * k_delta_avg) / ((N - 1) * (N - 2))

        self.links = [[] for _ in range(N)]
        self.triangles = [[] for _ in range(N)]

    def generate(self, seed=None):
        """Generates the topology using an optimized sampling approach."""
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

        self._generate_1_simplices()
        self._generate_2_simplices()
        return self.links, self.triangles

    def _generate_1_simplices(self):
        """Edge sampling from binomial."""
        print(f"Sampling edges with p_1 = {self.p_1:.8f}")
        num_edges = np.random.binomial(math.comb(self.N, 2), self.p_1)
        self.added_edges = set()

        node_list = range(self.N)
        while len(self.added_edges) < num_edges:
            i, j = random.sample(node_list, 2)
            edge = (i, j) if i < j else (j, i)

            if edge not in self.added_edges:
                self.added_edges.add(edge)
                self.links[i].append(j)
                self.links[j].append(i)

            print(f"\rEdges sampled: {len(self.added_edges)}/{num_edges}", end="")
        print()

    def _generate_2_simplices(self):
        """Triangle hyperedge sampling from binomial."""
        print(f"Sampling triangles with p_delta = {self.p_delta:.8f}")
        num_triangles = np.random.binomial(math.comb(self.N, 3), self.p_delta)
        added_triangles = set()

        node_list = range(self.N)
        while len(added_triangles) < num_triangles:
            nodes = tuple(sorted(random.sample(node_list, 3)))
            if nodes not in added_triangles:
                added_triangles.add(nodes)
                i, j, k = nodes
                self.triangles[i].append((j, k))
                self.triangles[j].append((i, k))
                self.triangles[k].append((i, j))

                for edge in [(i, j), (j, k), (i, k)]:
                    if edge not in self.added_edges:
                        self.added_edges.add(edge)
                        self.links[edge[0]].append(edge[1])
                        self.links[edge[1]].append(edge[0])

            print(f"\rTriangles sampled: {len(added_triangles)}/{num_triangles}", end="")
        print()
