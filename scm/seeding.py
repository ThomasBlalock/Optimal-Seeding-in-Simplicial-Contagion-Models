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
