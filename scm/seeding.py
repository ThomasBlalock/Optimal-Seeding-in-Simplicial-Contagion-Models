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
