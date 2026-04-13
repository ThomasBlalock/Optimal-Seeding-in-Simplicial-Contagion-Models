import numpy as np
import random


class SCMSimulator:
    """Synchronous-update SIS simulation over a simplicial complex.

    Each timestep: infected nodes recover with probability mu; susceptible
    nodes get infected with probability 1 - (1-beta)^m * (1-beta_delta)^n
    where m = infected neighbors via edges, n = infected triangles (both
    other members infected).
    """

    def __init__(self, links, triangles, initial_infected, beta, beta_delta, mu):
        self.links = links
        self.triangles = triangles
        self.N = len(links)

        self.beta = beta
        self.beta_delta = beta_delta
        self.mu = mu

        # 0 = Susceptible, 1 = Infected
        self.current_state = np.zeros(self.N, dtype=int)
        self.current_state[initial_infected] = 1

        self.rho_history = [np.mean(self.current_state)]

    def run(self, t_max):
        """Executes the simulation for t_max timesteps."""
        for t in range(t_max):
            self._step()
            if self.rho_history[-1] == 0.0:
                self.rho_history.extend([0.0] * (t_max - t - 1))
                break

        return self.rho_history

    def _step(self):
        next_state = np.copy(self.current_state)

        for i in range(self.N):
            if self.current_state[i] == 1:  # Infected
                if random.random() < self.mu:
                    next_state[i] = 0
            else:  # Susceptible
                m = self._count_infected_links(i)
                n = self._count_infected_triangles(i)

                p_inf = 1.0 - ((1.0 - self.beta) ** m) * ((1.0 - self.beta_delta) ** n)

                if random.random() < p_inf:
                    next_state[i] = 1

        self.current_state = next_state
        self.rho_history.append(np.mean(self.current_state))

    def _count_infected_links(self, i):
        return sum(self.current_state[neighbor] for neighbor in self.links[i])

    def _count_infected_triangles(self, i):
        count = 0
        for j, k in self.triangles[i]:
            if self.current_state[j] == 1 and self.current_state[k] == 1:
                count += 1
        return count
