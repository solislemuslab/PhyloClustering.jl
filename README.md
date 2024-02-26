# PhyloClustering.jl

[![Build Status](https://github.com/solislemuslab/PhyloClustering.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/solislemuslab/PhyloClustering.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://solislemuslab.github.io/PhyloClustering.jl/dev)
[![codecov.io](https://codecov.io/gh/YiboK/PhyloClustering.jl/branch/master/graph/badge.svg?token=AXGI6GHHCY)](http://codecov.io/gh/YiboK/PhyloClustering.jl)

## Overview

`PhyloClustering.jl` is a [Julia](http://julialang.org/) package to perform unsupervised learning on phylogenetic trees. The algorithms currently included are K-means, Hierarchical Clustering, Gaussian Mixture Model, and DBSCAN.

## Usage

`PhyloClustering.jl` is a Julia package, so the user needs to first install Julia, and then install the package.

To install Julia, follow the instructions in [here](https://julialang.org/downloads/).

To install the package, type inside Julia:
```julia
]
add PhyloClustering
```
## Help and errors

To get help, check the documentation [here](https://solislemuslab.github.io/PhyloClustering.jl/dev). Please report any bugs and errors by opening an
[issue](https://github.com/solislemuslab/PhyloClustering.jl/issues/new).

## Scripts

To check or get all the original data from tests, check the GitHub repository [here](https://github.com/YiboK/PhyloClustering-scripts).

## Citation

If you use `PhyloClustering.jl` in your work, we kindly ask that you cite the following paper: 
```
@article{kong_Tiley_solis-lemus_2023,
author = {Kong, Y., Tiley, G. P., and Sol'{i}s-Lemus, C.},
year = {2023},
title = {{Unsupervised learning of phylogenetic trees via split-weight embedding}},
url={https://arxiv.org/abs/2312.16074}
}
```

## License

`PhyloClustering.jl` is licensed under a
[MIT License](https://github.com/solislemuslab/PhyloClustering.jl/blob/master/LICENSE).

## Contributions

Users interested in expanding functionalities in `PhyloClustering.jl` are welcome to do so. See details on how to contribute in [CONTRIBUTING.md](https://github.com/solislemuslab/PhyloClustering.jl/blob/master/CONTRIBUTING.md).

