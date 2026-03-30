Optimal Seeding in Simplicial Contagion Models
By Thomas Blalock
March 2026

Introduction to the paper and why its an important work

Contagion refers to the process by which an endemic such as a virus, information, or an opinion spreads across a network. Traditional models of contagion rely on graphs where interactions are strictly dyadic, meaning that the contagion is modeled to traverse along the simple edges of the network. This approach is inadequate in the context of social contagion, where contagion through group interactions behaves fundamentally differently from contagion through pairwise interactions. For example, someone adopting a norm because other people independently convinced them to can not be equated to them adopting a norm because everyone in a group they participate in are all adopters. This is the difference between simple contagion and complex contagion.

To model this difference between pairwise interactions and group interactions, Iacopini 2019 formalizes different interactions as simplexes. Recalling the definition of a simplex, a k-simplex sigma is a set of k+1 vertices. A 0-simplex would be just a single node. A 1-simplex would be a pair of nodes with an edge connecting them. A 2-simplex or higher would be the set of k+1 nodes with a hyperedge connecting them in a group interaction. 

With this definition of group interactions as simplexes, Iacopini 2019 proposed the simplicial contagion model (SCM) as a modeling framework for contagion in social contexts, unifying simple and complex contagion into a single coherent modeling framework. The SCM treats different k-simplex interactions as containing fundamentally different processes by describing contagion as a stochastic process where the infection rates depend on the dimensionality of the interaction. In simpler terms, the infection rate of a 1-simplex (lambda) and that of a 2-simplex (lambda delta) are independent of each other. When modeled this way, group influence cannot be decomposed into independent pairwise edges.

Iacopini 2019 tested the SCM on several empirical and synthetic networks, yielding two primary conclusions: (i) higher-order interactions change the epidemic threshold transition from continuous to discontinuous and (ii) a bistable region emerges where both healthy and endemic states co-exist depending on the initial infection density. The SCM featuring a first-order phase transition explains well the “critical mass” observed to be required to trigger large-scale social change.

Overview of the reproduced result

To perform their experiments, they generated Erdos Renyi-style random graphs on a scale ranging from low rescaled infectivity to high rescaled infectivity, calculated as a function of the (insert variables) according to the formula: (insert formula) where (insert variable definitions). This rescaled infectivity score represents the potency of the contagion, meaning that a higher rescaled infectivity should correlate with a higher final infection density. They simulated the SCM over these synthetic networks at varying levels of 2-simplex infectivity rates (lambda delta). The results from these simulations display the discontinuous phase transition.

Given the use of Erdos Renyi-style random graph structures, they characterized the discontinuous phase transition using the mean field approach, which provided a theoretical convergence node density to check their simulations against. This yielded the following analytical solution for final node density: (insert mean field approach eq 1) where (insert variable definitions). This analytical solution correctly captures the phase transition since the final infection density is imaginary, and thus zero, when the rescaled infectivity is below the threshold but jumps to a high non-zero value above the rescaled infectivity threshold.

Figure 1a plots these simulations for varying values of rescaled infectivity and lambda delta, showing a discontinuous phase transition for each set of simulations using a non-zero lambda delta. But, when lambda delta is zero, meaning that the model is reduced to a simple contagion model without higher order interactions, the visualization displays a second order continuous phase transition.

Figure 1b plots simulations for varying initial densities of infected nodes over time with lambda fixed at 0.75 and lambda delta fixed at 2.5. For this configuration, the theoretical initial node density threshold is 0.2, plotted as a dashed line. These temporal evolutions of the infection density show the simulations starting above the threshold converging to the theoretical final infection density while the simulations with an initial infection density below the threshold converge to an infection density of zero. This shows the bistable region where the infection density is either zero or converges into a state where both infected and healthy nodes coexist.

(insert figures 1a and 1b)

Identify a criteria for success: how do you know when you’ve reproduced a result?

We aim to replicate their conclusion regarding the SCM yielding a discontinuous phase transition and the emergence of the bistable region dependent on initial node density. In order to perform this replication, we aim to recreate their simulations over synthetic 2-simplex networks, recreating figures 1a and 1b.

To recreate figure 1a, we must simulate the SCM for varying rescaled infectivities and values of lambda delta. Our simulation results must feature a discontinuous phase transition that approximately coincides with the mean field approach analytical solution as described in equation (insert equation number). Including the mean field approach solution validates the implementation of our random graph generation and simulation mechanics. Ensuring that our simulation results approximately match, not only the results in Iacopini 2019, but also the theoretical results as calculated by the mean field approach ensures that our models are not globally biased in either direction.

To recreate figure 1b, we must simulate the SCM over varying values of initial node density with a fixed rescaled infectivity and track the infection density over time. The results of these simulations must produce a bistable region, with simulations starting above the theoretical threshold as calculated by the mean field approach converging to the predicted final infection density and simulations below that threshold converging to zero. In these simulations, we expect the threshold to be sensitive to stochastic noise in the SCM relative to the size of the random network. This is because, on small networks, random chance may push the infection density of simulations starting near the threshold above the threshold. On larger networks, the noise will be less impactful to the overall network density, and thus the results will be less sensitive to noise.

Did it work? Try to understand why or why not.

It did work. Go over the intuition for why this happens and how you understand it. Talk about the effects of using an Erdos Renyi in that there probably is not much difference between random seeding and deliberate seeding. Explain the little bit of noise on the bistable graph and explain that it is just noise and that increasing the size of the network makes it cohere more to the theoretical results.

Our simulations, as described above, yielded approximately the same final infection densities as Iacopini 2019. Our simulations, using the same simulation and graph construction hyperparameters as Iacopini 2019, are shown in figures 2a and 2b. Figure 2a shows the discontinuous phase transition for each SCM simulation, and each simulation approximately matches the mean field approach analytical solution. Figure 2b shows the bistable region dependent on initial node density, with simulations starting above the density threshold converging to an endemic state and simulations starting below the density threshold converging to a healthy state. 

(insert figures 2a and 2b)

Notably, the simulation with an initial node density of (insert outlier node density), which is below the threshold, converged to the endemic state. One explanation is that this exception was due to random noise since the network was small (2,000 nodes). To test this explanation, we reran the simulation over larger random networks with 8,000 nodes. The results of these simulations, shown in figure 3, coincided with the theoretical results. Thus, we concluded that we had successfully replicated the results of Iacopini 2019.

(insert figure 3)

Reflect on the process, what did you learn by doing this task?

The process of replicating the results of Iacopini 2019 taught me to have redundant tests that my experiments are working correctly. I had an initial false start regarding the construction of the random graphs where I had neglected to add the simple edges that corresponded to the 2-simplex edges. My initial batch of simulations still featured the discontinuous phase transition. So, if that were my only validating check, then I might have missed the oversight in my code. However, visualising the results against the mean field approach analytical solution made the nature of the problem obvious. The final infection densities of my simulations were globally negatively biased. The redundant checks helped me catch a vital error. The observation also narrowed the location error substantially.

The act of building and running the experiments also made the limitations of their conclusions obvious, specifically regarding the use of Erdos-Renyi-style random graphs and random seeding. For me, a lot of the time, it is difficult to determine where there is opportunity for incremental improvements on existing methods. Instead, if a methodology does not exactly meet my needs, I tend to attempt to create my own methodologies from scratch. Recreating these experiments made it obvious to me how I could incrementally build upon their framework to solve problems relevant to my projects. I chose to replicate this paper because I figured it would be useful when building an algorithm for my app deciding to which users I would send notifications regarding a new event. Before replicating the results of this paper, I envisioned jumping straight to optimizing the placement of seed nodes, which would have caused me to skip steps in the research process. This replication made it more obvious to me that I should first seek to understand the properties of the SCM over networks that better resemble social networks than Erdos-Renyi random graphs, in my context with community structure and exponentially distributed degree. My intuition is that random seeding within communities or clusters biased by degree would dramatically lower the initial infection density threshold, providing a much better baseline for my eventual optimization algorithm. Forcing myself to actually replicate the paper rather than just read it gave me a more concrete direction for future research.




