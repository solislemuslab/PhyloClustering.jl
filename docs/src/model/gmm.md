# Gaussian Mixture Model (GMM)

[Gaussian Mixture Model (GMM)](https://scikit-learn.org/stable/modules/mixture.html#:~:text=A%20Gaussian%20mixture%20model%20is,Gaussian%20distributions%20with%20unknown%20parameters) a model-based probabilistic method that assumes data points are generated from a mixture of several Gaussian distributions with unknown parameters. It uses the Expectation Maximization
(EM) algorithm to update the parameters iteratively in order to optimize the log-likelihood of the data until
convergence.

```@docs
gmm_label
```

**Reference**

The implementation of *GMM* is provided by [`GaussianMixtures.jl`](https://github.com/davidavdav/GaussianMixtures.jl).