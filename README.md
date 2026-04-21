# Optimizing Seeding in Competitive Independent Capacity-Constrained Cascades

## Competitive Independent Capacity-Constrained Cascades (CIC3)

Many diffusion processes are only valuable up to a hard quota. For instance, social activities that leverage social contagion to recruit participants might only require a certain number of attendees before each subsequent attendee provides no value to the system or detracts value. A full social volleyball game only requires 12 people; adding more people to the game would detract from the net experience by causing some people to watch from the sidelines for periods of time. We refer to these types of contagions as **capacity-constrained cascades (C3s)**. Sometimes, many C3s are introduced simultaneously competing for attention with the goal of each one meeting its quota. We call this system **competitive independent capacity-constrained cascades (CIC3)**. Expanding on the previous example, if there were multiple social events occurring at the same time, then each event would be competing for attendees in the CIC3 process.

## Global CIC3 Evaluation - Deriving Time-Discounted Global Attainment

When evaluating CIC3 processes, we want to capture a few aspects of the system. The first is global attainment: we want to know how much of each contagion's quota was filled. The following definition matches this goal by only including infected nodes towards the numerator up to a given contagion's quota. 

Let $\mathcal{C}$ be the set of contagions and let $C_i \in \mathcal{C}$ denote contagion $i$. Let $I_i^{\text{raw}}$ be the set of nodes that were infected by $C_i$, let $t(v)$ denote the infection timestep of node $v \in I_i^{\text{raw}}$, and let $Q_i$ be the quota for $C_i$. Define the capped infection count:

$$K_i = \min(|I_i^{\text{raw}}|, Q_i)$$

Then per-contagion attainment ($A_i$) and global attainment ($A_g$) are:

$$A_i = \frac{K_i}{Q_i} \in [0,1]$$

$$A_g = \frac{1}{|\mathcal{C}|} \sum_{C_i \in \mathcal{C}} A_i \in [0,1]$$

This formulation assumes that the values of each contagion meeting its quota are equal. If they were not equal, it would be trivial to redefine global attainment as a weighted average rather than the mean.

In most use cases, the speed of contagion would also be of interest. Infecting a node at an earlier timestep is more valuable than infecting a node at a later timestep. Assuming we have some value function $V(t)$ that maps an infection timestep to a value in $[0,1]$, we can reformulate our metric to include a time discount.

Let $t_{i,(1)} \le t_{i,(2)} \le \dots \le t_{i,(|I_i^{\text{raw}}|)}$ be the sorted infection times for nodes in $I_i^{\text{raw}}$. Define the time-discounted capped count as the sum of values for the earliest $\min(Q_i,|I_i^{\text{raw}}|)$ infections, truncated at $Q_i$:

$$K_i^{\text{td}} = \min \left( \sum_{k=1}^{\min(Q_i, |I_i^{\text{raw}}|)} V(t_{i,(k)}), Q_i \right)$$

Then:

$$A_i^{\text{td}} = \frac{K_i^{\text{td}}}{Q_i} \in [0,1]$$

$$A_g^{\text{td}} = \frac{1}{|\mathcal{C}|} \sum_{C_i \in \mathcal{C}} A_i^{\text{td}} \in [0,1]$$

## Comparison with Prior Metrics

Prior metrics for network contagion do not capture both a hard quota and a time-discounted value under that quota.

| Metric | What it measures | Why it is inadequate for CIC3 |
| :--- | :--- | :--- |
| **Final epidemic size / prevalence** | Total nodes ever infected | Monotone in infections. No hard quota $Q_i$ where value plateaus; counts infections beyond $Q_i$ as benefit. No time weighting. |
| **Basic reproduction number $R_0$** | Expected secondary infections per primary | No notion of a quota. No time-value function $V(t)$. |
| **Expected spread $\sigma(S)$** | Expected activated nodes from seed set $S$ | Objective increases with every additional adoption past $Q_i$; does not enforce $\min(|I_i^{\text{raw}}|, Q_i)$. |
| **Structural virality** | Average distance in diffusion tree | Describes shape of spread, not whether $Q_i$ was met or how early infections occurred. |
| **Time to $k$ infections** | Speed to reach a fixed count $k$ | Measures latency but does not cap value at $Q_i$ and does not apply a graded $V(t)$ to all infections within the quota. |

## Experiments

### CIC3 Seeding Strategy Comparison

The `experiments/cic3_seeding_strategies.ipynb` notebook compares five seeding strategies across four network topologies under the CIC3 evaluation metric. Each (topology, strategy) cell runs `NUM_TRIALS` independent simulations and aggregates the results.

**Strategies.** Three baselines plus two new clustering-based strategies:

- **Random** — uniformly sampled disjoint seed sets per contagion.
- **High Degree** — top-degree nodes assigned round-robin across contagions.
- **High 2-Simplex** — top-triangle-participation nodes assigned round-robin across contagions.
- **Louvain** — runs a single-level Louvain community detection pass on the edge list, ranks communities by internal edge-endpoint density, and fills each of the top $C$ communities with the highest-internal-degree nodes in that community. Falls back to Farthest-First if Louvain finds fewer than $C$ communities.
- **Farthest-First** — picks $C$ BFS-farthest centers (explicitly avoiding high-degree hub nodes as centers) and grows a BFS-ball of `num_seeds_per_contagion` around each center, preferring low-degree community members within each BFS layer so hubs are not pulled in.

Both clustering-based strategies discover community structure directly from the graph with no external labels — they use the same `(N, num_seeds_per_contagion, links, triangles)` constructor as the baselines.

**Topologies.** `RSC` (random simplicial complex, Erdős–Rényi-style), `BA` (Barabási–Albert with a simplicial extension), `SBM` (a small stochastic block model), and `SBM-Manual` — an SBM configured to stress clustering-based strategies:

- $N = 500$, $K = 53$ communities with sizes `[9]*52 + [32]`.
- 52 small communities are full cliques (`p_intra = 1.0`), with very low cross-community bleed (`p_inter = 0.003`).
- One hub community (size 32) connects densely to every other community (`p_inter = 0.15`) and has moderate internal density (`p_intra = 0.5`). Hub nodes dominate the degree distribution.
- Triangle probabilities: `[0.2]*52 + [0.04]`. Realized $\langle k \rangle \approx 19.05$, $\langle k_\Delta \rangle \approx 6.54$.

**Metric.** Global attainment $A_g$ and time-discounted global attainment $A_g^{\text{td}}$ with $V(t) = e^{-0.05 t}$. Quotas are set so $\sum_i Q_i = N$ (no slack), with $C = 50$ contagions and $Q_i = N/C = 10$.

### Charts

**`figures/cic3_strategy_comparison.png`** — Grouped bar chart comparing the five strategies across the four topologies. Solid bars show mean $A_g^{\text{td}}$ across trials; faint bars behind show mean $A_g$. Error bars are trial-to-trial standard deviation.

![CIC3 strategy comparison](figures/cic3_strategy_comparison.png)

**`figures/cic3_per_contagion_distribution.png`** — Per-topology boxplots of the pooled per-contagion $A_i^{\text{td}}$ values across all trials and contagions. One boxplot group per topology, with one box per strategy showing the median, IQR, whiskers, and outliers of individual contagion attainments.

![CIC3 per-contagion distribution](figures/cic3_per_contagion_distribution.png)

**`figures/cic3_per_contagion_pdf.png`** — Per-topology overlapping probability density curves of per-contagion $A_i^{\text{td}}$ built from `np.histogram(density=True)` (one curve per strategy per topology). Designed to surface shape features that boxplots flatten — bimodality, long lower tails, and saturation spikes at $A_i^{\text{td}} = 1$.

![CIC3 per-contagion PDF](figures/cic3_per_contagion_pdf.png)

## Notes

The high-degree seeding strategies tended to perform the best, except for over one topology: the twitter mutual network. Over that network, random seeindg performed the best. The community based seeding and furthest distance-based strategy did not work that well on any topology. These results were very odd, especially the random performing the best on the twitter mutual network. I statistically showed that the performance of random was statistically significantly better than the other strategies. Then, I examined various distributions in the graph topologies to try to find what made the twitter mutual network special. The primary difference between the twitter mutual network and the other topologies was that it was a lot sparser. It may be that degree-centrality seeding is the optimal strategy when the network is too dense for community structure to adequately contain contagion.

Experiment idea. Maybe we take the SBM and define some variables to tweak:
- Inter Community Connection: varying from isolated communities that only connect to hubs to no community structure at all.
- Intra Community Connection: varying from fully connected communities to sparsely connected communities
- Hubbyness: Fix this. We're investigating community structure and sparseness, not degree distribution
With this, make a heatmap with inter community conenction on the y axis and intra on the x axist and the performance of our selected seedin strategies as the heat value. Do this for each seeding strategy. Could also fix one or the other, sweer the non-fixed one, and plot all the seeding strategy's performance values on the same plot.

Update: The experiments did not show the expected pattern with random beating our high degree. Maybe we need to test a topology with a power law degree distribution and community structure?

Update: We tried the tests on a popularity-similarity topology with community structure and we were still unable to replicate the anomoly noticed on the twitter mutual network.

Update: We tried doing it on a network with super super sparse connections between communities and it worked finally, showing random to work better. The explination is obvious: using high-degree seeding makes certain parts of the network unreachable because they're disconnected or maybe blocked by a different contagion over a bridge.