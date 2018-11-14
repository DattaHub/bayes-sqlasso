# Adapting to Heavy Tails and Sparsity: Bayesian Square-root Lasso and Dirichlet-Sqrt-Lasso

Authors: M. Abba, J. Datta. 

## Abstract 
Abstract. The Lasso and the global-local shrinkage priors, gold-standards in the frequentist and Bayesian paradigms, critically depend on learning the error variance. This causes a lack of scale invariance and adaptability to heavy-tailed data. The Square-root Lasso (Belloni et al., 2011) attempts to correct this by using the $\ell_1$ norm on both the likelihood and the penalty for the objective function. In contrast, there is essentially no method for uncertainty quantification or automatic parameter tuning via a formal Bayesian treatment of an unknown error distribution. Here wpe propose a formal Bayesian solution to these problems, called the Bayesian Lasso and its extension Dirichlet-Square-Root-Lasso that achieve scale invariance and robpustness to heavy tails while maintaining computational efficiency. The Bayes Lasso leads to uncertainty quantification by yielding standard error estimates and credible sets for the underlying parameters. Furthermore, the hierarchical model leads to an automatic tuning of the penalty parameter using a full Bayes or empirical Bayes approach, avoiding any ad-hoc choice over a grid. An astonishing outcome is that these methods adapt to the entire range of sparsity, i.e handle ultra-sparse as well dense parameters, unlike any other global-local priors. We provide an efficient Gibbs sampling scheme based on Gaussian scale mixture identities. Performance on real and simulated data exhibit excellent small sample
properties and we establish some theoretical guarantees.
