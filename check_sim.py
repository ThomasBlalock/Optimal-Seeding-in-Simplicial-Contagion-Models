import sys
from scm.generators import TwitterEgoGenerator
from scm.seeding import RandomSeeding, HighDegreeSeeding, MultiRandomSeeding, MultiHighDegreeSeeding
from scm.cic3_simulator import CIC3Simulator
from scm.analysis import attainment

print("Generating Twitter Mutual...")
gen = TwitterEgoGenerator(p_promote=0.01, edge_mode="mutual", data_dir="Optimal-Seeding-in-Simplicial-Contagion-Models/data")
gen.generate()

N = gen.N
C = 50
Q = N // C
quotas = [Q] * C
num_seeds = [max(1, Q // 5)] * C

print(f"N={N}, C={C}, Quota={Q}, Seeds per contagion={num_seeds[0]}")

print("Running Random Seeding...")
strat_random = MultiRandomSeeding(N, num_seeds, gen.links, gen.triangles)
seeds_random = strat_random.seed()
sim_random = CIC3Simulator(gen.links, gen.triangles, seeds_random, [0.04]*C, [0.03]*C, quotas)
sim_random.run(t_max=200)
uninfected_random = sim_random.infected_by.tolist().count(-1)
A_i_r, A_g_r = attainment(sim_random.infected_by, quotas)
print(f"Random: {uninfected_random} uninfected nodes. Global Attainment: {A_g_r:.4f}")

print("Running High Degree Seeding...")
strat_hd = MultiHighDegreeSeeding(N, num_seeds, gen.links, gen.triangles)
seeds_hd = strat_hd.seed()
sim_hd = CIC3Simulator(gen.links, gen.triangles, seeds_hd, [0.04]*C, [0.03]*C, quotas)
sim_hd.run(t_max=200)
uninfected_hd = sim_hd.infected_by.tolist().count(-1)
A_i_hd, A_g_hd = attainment(sim_hd.infected_by, quotas)
print(f"High Degree: {uninfected_hd} uninfected nodes. Global Attainment: {A_g_hd:.4f}")

