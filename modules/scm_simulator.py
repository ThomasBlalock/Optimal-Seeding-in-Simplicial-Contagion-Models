import numpy as np
import random


class RandomSeeding:
    def __init__(self, N, rho_0):
        self.N = N
        self.rho_0 = rho_0

    def seed(self):
        num_initial_infected = int(self.N * self.rho_0)
        initial_infected = random.sample(range(self.N), num_initial_infected)
        return initial_infected


class SCMSimulator:
    def __init__(self, links, triangles, initial_infected, beta, beta_delta, mu):
        """Initializes the simulator with the given parameters.
        Parameters:
        - links: List of lists, where links[i] contains the neighbors of node i (1-simplices).
        - triangles: List of lists, where triangles[i] contains tuples (j, k) representing 2-simplices that include node i.
        - beta: Infection probability for 1-simplices.
        - beta_delta: Infection probability for 2-simplices.
        - mu: Recovery probability.
        - initial_infected: List of node indices that are initially infected.
        """
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
            if self.rho_history[-1] == 0.0: # Absorbing state
                # Pad for output len consistency & break early
                self.rho_history.extend([0.0] * (t_max - t - 1))
                break
                
        return self.rho_history

    def _step(self):
        next_state = np.copy(self.current_state)
        
        for i in range(self.N):
            if self.current_state[i] == 1: # Infected
                if random.random() < self.mu:
                    next_state[i] = 0
            else: # Susceptible
                m = self._count_infected_links(i)
                n = self._count_infected_triangles(i)
                
                # Total infection probability
                p_inf = 1.0 - ((1.0 - self.beta)**m)*((1.0 - self.beta_delta)**n)
                
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