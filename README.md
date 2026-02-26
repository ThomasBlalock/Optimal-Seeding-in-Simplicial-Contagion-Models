# Replication Plan: Simplicial Contagion Model (SCM)

## Phase 1: Substrate Generation (The Network Topology)
To replicate the findings, you first need a controllable environment. [cite_start]The authors use a Random Simplicial Complex (RSC) model[cite: 203]. 

**1. Build the RSC Generator Module**
Create a standalone class or module that generates an RSC. This isolates the topology from the contagion dynamics.
* [cite_start]**Target Parameters:** Set the number of nodes $N=2000$, average link degree $\langle k \rangle \simeq 20$, and average 2-simplex degree $\langle k_\Delta \rangle \simeq 6$[cite: 213, 282].
* [cite_start]**Edge Creation:** Connect any pair of vertices $(i, j)$ with probability $p_1$[cite: 205].
* [cite_start]**Triangle Creation:** Connect any triplet of vertices $(i, j, k)$ with probability $p_\Delta$[cite: 206].
* **Parameter Tuning Equations:** To achieve the target degrees, calculate your probabilities using the exact equations from the text:
    * [cite_start]$p_1 = \frac{\langle k \rangle - 2\langle k_\Delta \rangle}{N - 1 - 2\langle k_\Delta \rangle}$[cite: 403].
    * [cite_start]$p_\Delta = \frac{2\langle k_\Delta \rangle}{(N-1)(N-2)}$[cite: 409].
* **Data Structures:** Store the network as both an adjacency list (for fast 1-simplex neighbor lookups) and a list of triplets/sets (for fast 2-simplex lookups). 

## Phase 2: The SCM Simulator (The Dynamics)
Design the contagion simulation as a state machine that takes a graph object (from Phase 1) and an initial seed configuration as inputs. This modularity is what will allow you to inject your influence maximization algorithms later.

**1. Define the State and Parameters**
* [cite_start]**Node State:** Track a binary state variable $x_i(t) \in \{0,1\}$ for each node, where 0 is Susceptible and 1 is Infectious[cite: 91, 92].
* [cite_start]**Global Parameters:** * $\mu$: Node-independent recovery probability[cite: 100, 101].
    * [cite_start]$\beta$: Probability per unit time of infection via a 1-simplex (link)[cite: 95, 96].
    * [cite_start]$\beta_\Delta$: Probability per unit time of infection via a 2-simplex (full triangle)[cite: 96].
    * [cite_start]$\lambda = \beta\langle k \rangle / \mu$: Rescaled link infectivity[cite: 113].
    * [cite_start]$\lambda_\Delta = \beta_\Delta\langle k_\Delta \rangle / \mu$: Rescaled triangle infectivity[cite: 113].

**2. Implement the Timestep Update Logic**
At each discrete timestep $t$, update the states based on the following transition probabilities:
* [cite_start]**Recovery ($I \rightarrow S$):** Infected nodes recover with probability $\mu$[cite: 88, 100, 101].
* [cite_start]**Link Infection ($S \rightarrow I$):** A susceptible node gets infected with probability $\beta$ at each timestep through each link connected to an infected neighbor[cite: 85].
* [cite_start]**Triangle Infection ($S \rightarrow I$):** A susceptible node gets infected with probability $\beta_\Delta$ from a 2-face if the two other nodes in the 2-simplex are already infected[cite: 87, 96]. 
*(Note: Be careful to calculate the aggregate probability of infection for a susceptible node at time $t$ by evaluating all its active 1-simplex and 2-simplex channels before applying the state change, avoiding order-of-operation bias).*

**3. Implement the Observers**
* [cite_start]**Macroscopic Order Parameter:** Track the density of infectious nodes over time: $\rho(t) = \frac{1}{N}\sum_{i=1}^N x_i(t)$[cite: 94].
* [cite_start]**Stationary Density ($\rho^*$):** Build a function to average $\rho(t)$ over the last 100 time-steps after the system reaches an absorbing state or a stationary state[cite: 109].

## Phase 3: Analytical Baseline
Before running simulations, implement the Mean-Field (MF) equations. Plotting these first gives you a definitive target to unit-test your simulator against.

**1. Calculate Steady States**
[cite_start]Implement the analytical solution for the stationary density (Eq. 4)[cite: 290, 292]:
[cite_start]$$\rho_{2\pm}^* = \frac{\lambda_\Delta - \lambda \pm \sqrt{(\lambda-\lambda_\Delta)^2 - 4\lambda_\Delta(1-\lambda)}}{2\lambda_\Delta}$$ [cite: 293]
* [cite_start]Plot this for varying $\lambda$ and $\lambda_\Delta$ to recreate the red theoretical curves seen in Fig 3a[cite: 283, 285].

## Phase 4: Replication Execution
Now, plug the components together to recreate the specific figures. 

**1. Replicating Figure 3a (Phase Transition)**
* [cite_start]**Setup:** Generate new RSC instances for different runs[cite: 213]. Fix $\mu$ and sweep $\lambda \in [0, 2]$. 
* [cite_start]**Test Cases:** Run simulations for $\lambda_\Delta \in \{0, 0.8, 2.5\}$[cite: 214]. 
* [cite_start]**Seeding:** Start with a randomly placed initial density $\rho_0$[cite: 110, 213]. 
* [cite_start]**Output:** Plot the average stationary density $\rho^*$ against $\lambda$ to demonstrate the continuous transition at $\lambda_\Delta = 0.8$ and the discontinuous transition (with hysteresis) at $\lambda_\Delta = 2.5$[cite: 217, 218]. 

**2. Replicating Figure 3b (Bistability and Critical Mass)**
* [cite_start]**Setup:** Fix the parameters within the bistable region: $\lambda = 0.75$ and $\lambda_\Delta = 2.5$[cite: 288].
* [cite_start]**Sweep Initial Density:** Vary the initial seed density $\rho_0 \equiv \rho(0)$ using different random subsets of nodes[cite: 289].
* [cite_start]**Output:** Plot the temporal evolution of $\rho(t)$ up to $t=500$[cite: 254, 274, 279, 280]. 
* [cite_start]**Validation:** Verify that trajectories starting below the MF unstable branch $\rho_{2-}^*$ collapse to 0, and those starting above it rise to the endemic state[cite: 290, 310, 314]. 

## Phase 5: Groundwork for Extension
Because Phase 2 decoupled the seeding from the dynamics, moving to your PhD extension requires no changes to the SCM engine. 
* You simply replace the `random_seed(rho_0)` function with your custom algorithmic seeders (e.g., `greedy_simplex_centrality_seed(rho_0)`).
* You can then run the exact same Fig 3b analysis to see if your optimal seeding effectively lowers the required critical mass ($\rho_{2-}^*$) compared to random seeding.