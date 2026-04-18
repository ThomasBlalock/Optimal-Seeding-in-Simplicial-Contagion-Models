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