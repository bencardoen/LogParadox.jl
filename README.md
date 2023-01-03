# LogParadox

A project to illustrate how you can obtain paradoxical pattern inversions when applying hypothesis tests, conditional on a log transform of your data.

[![Coverage](https://codecov.io/gh/bencardoen/LogParadox.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/bencardoen/LogParadox.jl)

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/bencardoen/LogParadox.jl/tree/main.svg?style=svg&circle-token=304e0f4d40f0fdb0363572f8fabf8ee73334ebfd)](https://dl.circleci.com/status-badge/redirect/gh/bencardoen/LogParadox.jl/tree/main)

## Example
Suppose you have a vector of values X (box plot in blue).
This example shows we can find X' such that E[X] > E[X'], yet E[log(X)] <  E[log(X')].
In other words, there exists distributions, especially in long tail data, where small differences between the 2 datasets can induce a 'X is greater than Y' conclusion, yet in log scale, report a 'X is smaller than Y' conclusion.
We show that you can get this effect with as little as 5% of the data differing.
The necessary and sufficient conditions are derived in the paper, but our API allows you to test to see if you data is vulnerable or not.
![example](figures/interactivelp.gif)

### Effect on significance testing
Using a non-parametric hypothesis test, at 5% of data modified you can induce a strong effect consistly achieving significance.
Note that this does not try to reinforce the flawed idea that significance in isolation is a valid finding, rather serve as a cautionary tale that inducing significant inversions is fairly easily done.
![example](figures/pvals50_200.png)

## Installation
- Get [Julia](https://julialang.org/learning/getting-started/)
```bash
julia
julia 1.x>using Pkg; Pkg.add(url="https://github.com/bencardoen/ERGO.jl.git")
julia 1.x>using Pkg; Pkg.add(url="https://github.com/bencardoen/SPECHT.jl.git")
julia 1.x>using Pkg; Pkg.add(url="https://github.com/bencardoen/LogParadox.jl.git")
julia 1.x>Pkg.test("DataCurator")
julia 1.x>using DataCurator
```
or
```bash
git clone https://github.com/bencardoen/LogParadox.jl.git
cd LogParadox.jl
julia
julia 1.x>using Pkg; Pkg.activate(); Pkg.resolve(); Pkg.instantiate(); Pkg.test()
julia 1.x>using DataCurator
```

## Usage

See [scripts](scripts/gif.jl) for example illustrations that use the API.

## Cite
```bibtex
TODO
```
