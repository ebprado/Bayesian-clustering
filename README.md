# Bayesian clustering
In order to have a deeper understanding about some basic Bayesian clustering methods, such as Finite Mixture of Normals (FMN) and Product Partition Models (PPM), I decided to implement them in R as simple as possible so that I could grasp.

## Finite Mixture of Normals
Suppose we have a set of samples X_{1}, X_{2}, ..., X_{n$ can be modelled as
\begin{align}
    p(X_{i}, | \mu_{1:k}, \tau_{1:k}, q_{1:k}) = \sum_{j}^{k} q_{j} \mbox{N}(\mu_{j}, \tau_{j}^{-1}),
\end{align}
where \mbox{N}(.) denotes the Normal distribution and $q_{j}$ represents the weight of the j-th component, with $q_{j} > 0$ and $\sum_{j=1}^{k} q_{j} = 1$. In addition, let's introduce the latent variable $Z_{1:n}$ to induce the mixture. Thus, we have that $X_{i}|Z_{i} = j \sim \mbox{N}(\mu_{j}, \tau_{j})$. Given the introduction of the latent variable, we can rewrite the likelihood function as
