import numpy as np


def mf_rho(lam, lam_d):
    """Mean-field analytical steady state of the SCM.

    Returns the real positive root of:
        rho = (lam_d - lam +/- sqrt((lam - lam_d)^2 - 4*lam_d*(1 - lam))) / (2*lam_d)
    or 0 when the discriminant is negative (below threshold).
    """
    if lam_d == 0:
        return 1 - 1 / lam

    discriminant = (lam - lam_d) ** 2 - 4 * lam_d * (1 - lam)
    sqrt_term = np.sqrt(discriminant.astype(complex))
    rho = (lam_d - lam + sqrt_term) / (2 * lam_d)
    rho = np.where(np.imag(rho) != 0, 0, rho)

    return rho.real


def steady_state_rho(rho_history, t_avg=100):
    """Compute steady-state infection density from a simulation history."""
    if rho_history[-1] == 0.0:
        return 0.0
    return np.mean(rho_history[-t_avg:])


def attainment(infected_by, quotas):
    """CIC3 global attainment (no time discount).

    infected_by: array of shape (N,) with values in {-1, 0..C-1}
        (-1 = never infected, otherwise the contagion index).
    quotas: array-like of length C (positive ints).

    Returns (A_i, A_g):
        A_i: np.ndarray of shape (C,), per-contagion attainment in [0,1]
        A_g: float, mean over contagions
    """
    infected_by = np.asarray(infected_by)
    quotas = np.asarray(quotas, dtype=float)
    C = len(quotas)

    A_i = np.zeros(C)
    for c in range(C):
        raw = int((infected_by == c).sum())
        A_i[c] = min(raw, quotas[c]) / quotas[c]

    return A_i, float(A_i.mean())


def time_discounted_attainment(infected_by, infection_times, quotas, V):
    """CIC3 time-discounted global attainment.

    Sums V(t) over the earliest min(Q_i, |I_i^raw|) infections per
    contagion, truncates at Q_i, divides by Q_i.

    V: callable int -> float in [0,1] (value of infecting a node at step t).

    Returns (A_i_td, A_g_td) with the same shapes as `attainment`.
    """
    infected_by = np.asarray(infected_by)
    infection_times = np.asarray(infection_times)
    quotas = np.asarray(quotas, dtype=float)
    C = len(quotas)

    A_i_td = np.zeros(C)
    for c in range(C):
        mask = infected_by == c
        times = np.sort(infection_times[mask])
        Q_i = int(quotas[c])
        k = min(Q_i, len(times))
        if k == 0:
            A_i_td[c] = 0.0
            continue
        values = np.array([V(int(t)) for t in times[:k]])
        K_i_td = min(values.sum(), quotas[c])
        A_i_td[c] = K_i_td / quotas[c]

    return A_i_td, float(A_i_td.mean())


def exponential_decay(rate):
    """V(t) = exp(-rate * t)."""
    def V(t):
        return float(np.exp(-rate * t))
    return V


def linear_decay(T):
    """V(t) = max(0, 1 - t/T)."""
    def V(t):
        return max(0.0, 1.0 - t / T)
    return V


def no_decay():
    """V(t) = 1 for all t."""
    def V(t):
        return 1.0
    return V


def deadweight(infected_by, quotas):
    """CIC3 deadweight: nodes infected beyond the quota per contagion.

    Deadweight for contagion c = max(0, |I_c| - Q_c) where |I_c| is the
    raw count of nodes infected by c and Q_c is its quota.

    Args:
        infected_by: array of shape (N,) with values in {-1, 0..C-1}.
        quotas: array-like of length C (positive ints).

    Returns:
        (D_i, D_g):
            D_i: np.ndarray of shape (C,), per-contagion deadweight count
            D_g: float, mean deadweight across contagions
    """
    infected_by = np.asarray(infected_by)
    quotas = np.asarray(quotas, dtype=int)
    C = len(quotas)

    D_i = np.zeros(C, dtype=int)
    for c in range(C):
        raw = int((infected_by == c).sum())
        D_i[c] = max(0, raw - quotas[c])

    return D_i, float(D_i.mean())


def penetration(links, infected_by, seed_sets):
    """CIC3 penetration: average BFS hop distance from seed nodes
    through the subgraph induced by each contagion's infected nodes.

    For contagion c, BFS starts from all seed nodes of c and expands
    only through edges where both endpoints are infected by c.
    Penetration P_c = mean hop distance over all non-seed nodes
    infected by c. If contagion c infected only seed nodes (or none),
    P_c = 0.

    Args:
        links: adjacency list, links[i] = list of neighbor indices.
        infected_by: array of shape (N,) with values in {-1, 0..C-1}.
        seed_sets: list of C lists, seed_sets[c] = list of seed node
            indices for contagion c.

    Returns:
        (P_i, P_g):
            P_i: np.ndarray of shape (C,), per-contagion penetration
                (mean hop distance from seeds). 0 if no non-seed
                infections.
            P_g: float, mean penetration across contagions
    """
    from collections import deque

    infected_by = np.asarray(infected_by)
    C = len(seed_sets)

    P_i = np.zeros(C)
    for c in range(C):
        seeds = seed_sets[c]
        if not seeds:
            continue

        # BFS restricted to nodes infected by contagion c
        dist = {}
        queue = deque()
        for s in seeds:
            dist[s] = 0
            queue.append(s)

        while queue:
            node = queue.popleft()
            d = dist[node]
            for nb in links[node]:
                if nb not in dist and infected_by[nb] == c:
                    dist[nb] = d + 1
                    queue.append(nb)

        # Average over non-seed infected nodes
        non_seed_hops = [dist[n] for n in dist if dist[n] > 0]
        if non_seed_hops:
            P_i[c] = float(np.mean(non_seed_hops))

    return P_i, float(P_i.mean())
