import numpy as np
import random


class CIC3Simulator:
    """Synchronous-update SI simulation of competing capacity-constrained
    contagions (CIC3) on a simplicial complex.

    Semantics:
    - Exclusive states: a node can host at most one contagion. Once
      infected by contagion c, it is locked out of all others.
    - SI dynamics: no recovery.
    - Quotas cap only the evaluation metric; contagions keep spreading
      after their quota is met. Simulation terminates when no
      susceptibles remain or at t_max. If stop_on_all_quotas_met is True,
      also stops when all quotas are met (legacy behavior).

    State encoding: current_state[i] is 0 for susceptible or c+1 for a node
    infected by contagion c (1-indexed so the default 0 means susceptible).
    infected_by[i] is -1 or c in [0, C-1] for the metric.

    Per-step transition for each susceptible i:
      1. For each contagion c compute p_c = 1 - (1-beta_c)^m_c * (1-beta_delta_c)^n_c
         where m_c = edge-neighbors with state == c+1 and n_c = triangles
         where both other members have state == c+1.
      2. Independently roll each contagion; let S = {c : roll succeeded}.
      3. |S|=0: remain susceptible. |S|=1: infect with that contagion.
         |S|>1: break ties via tie_break.
    """

    def __init__(
        self,
        links,
        triangles,
        initial_infected_per_contagion,
        betas,
        beta_deltas,
        quotas,
        tie_break="uniform",
        stop_on_all_quotas_met=False,
    ):
        self.links = links
        self.triangles = triangles
        self.N = len(links)

        self.betas = np.asarray(betas, dtype=float)
        self.beta_deltas = np.asarray(beta_deltas, dtype=float)
        self.quotas = np.asarray(quotas, dtype=int)
        self.C = len(self.betas)

        if not (len(self.beta_deltas) == self.C == len(self.quotas)):
            raise ValueError(
                "betas, beta_deltas, and quotas must have the same length"
            )
        if len(initial_infected_per_contagion) != self.C:
            raise ValueError(
                f"initial_infected_per_contagion has {len(initial_infected_per_contagion)} "
                f"sets but there are {self.C} contagions"
            )
        if tie_break not in ("uniform", "intensity"):
            raise ValueError(
                f"unknown tie_break = {tie_break!r}; expected 'uniform' or 'intensity'"
            )
        self.tie_break = tie_break
        self.stop_on_all_quotas_met = stop_on_all_quotas_met

        self.current_state = np.zeros(self.N, dtype=int)
        self.infected_by = np.full(self.N, -1, dtype=int)
        self.infection_times = np.full(self.N, -1, dtype=int)

        seen = set()
        for c, seeds in enumerate(initial_infected_per_contagion):
            dup = seen.intersection(seeds)
            if dup:
                raise ValueError(
                    f"Seed sets not disjoint: contagion {c} shares nodes "
                    f"{sorted(dup)} with an earlier contagion"
                )
            seen.update(seeds)
            for i in seeds:
                self.current_state[i] = c + 1
                self.infected_by[i] = c
                self.infection_times[i] = 0

        self.rho_history = [self._per_contagion_rho()]

    def run(self, t_max):
        """Execute up to t_max timesteps. Returns rho_history as a
        (T, C) array, padded to length t_max + 1 if an early-stop
        condition fires."""
        for t in range(1, t_max + 1):
            self._step(t)

            if self._no_susceptibles() or (
                self.stop_on_all_quotas_met and self._all_quotas_met()
            ):
                final = self.rho_history[-1]
                while len(self.rho_history) < t_max + 1:
                    self.rho_history.append(final)
                break

        return np.array(self.rho_history)

    def _step(self, t):
        next_state = np.copy(self.current_state)
        next_infected_by = np.copy(self.infected_by)
        next_infection_times = np.copy(self.infection_times)

        for i in range(self.N):
            if self.current_state[i] != 0:
                continue

            p = self._per_contagion_infection_prob(i)
            winners = [c for c in range(self.C) if random.random() < p[c]]

            if not winners:
                continue
            if len(winners) == 1:
                c = winners[0]
            else:
                c = self._break_tie(winners, p)

            next_state[i] = c + 1
            next_infected_by[i] = c
            next_infection_times[i] = t

        self.current_state = next_state
        self.infected_by = next_infected_by
        self.infection_times = next_infection_times
        self.rho_history.append(self._per_contagion_rho())

    def _per_contagion_infection_prob(self, i):
        p = np.zeros(self.C)
        for c in range(self.C):
            m = self._count_infected_links(i, c)
            n = self._count_infected_triangles(i, c)
            p[c] = 1.0 - ((1.0 - self.betas[c]) ** m) * (
                (1.0 - self.beta_deltas[c]) ** n
            )
        return p

    def _count_infected_links(self, i, c):
        label = c + 1
        return sum(1 for j in self.links[i] if self.current_state[j] == label)

    def _count_infected_triangles(self, i, c):
        label = c + 1
        count = 0
        for j, k in self.triangles[i]:
            if self.current_state[j] == label and self.current_state[k] == label:
                count += 1
        return count

    def _break_tie(self, winners, p):
        if self.tie_break == "uniform":
            return random.choice(winners)
        # intensity: weight by -log(1 - p_c), the underlying Poisson rate
        weights = np.array([-np.log(1.0 - p[c]) for c in winners])
        weights = weights / weights.sum()
        return int(np.random.choice(winners, p=weights))

    def _per_contagion_rho(self):
        rho = np.zeros(self.C)
        for c in range(self.C):
            rho[c] = (self.infected_by == c).sum() / self.N
        return rho

    def _all_quotas_met(self):
        for c in range(self.C):
            if (self.infected_by == c).sum() < self.quotas[c]:
                return False
        return True

    def _no_susceptibles(self):
        return not (self.current_state == 0).any()
