# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# Copyright 2022-3, Ben Cardoen
using Plots
using LogParadox
using Distributions
default(legend = false)
using Random
using StatsPlots
default(legend = false)
Random.seed!(42)

### A small script to plot a gif showing how the paradox can be induced on Log-Normal data

### Set parameters
μ, σ, N, M = 7, 2, 1000, 50

### Generate in silico data
xs = exp.(rand(Normal(μ, σ), N))
sx = [sample(xs, M) for _ in 1:M]

### Transform
warpedx = minmaxreplace.(sx)

### Compute means
gs = gm.(sx)
ms = am.(sx)

gsp = gm.(warpedx)
msp = am.(warpedx)

### Plot
p=boxplot([gs,gsp])
q=boxplot([ms,msp])
Plots.plot(p, q, xticks=[])

### Same thing, but with a gif to illustrate change
S, T = 200, 100
sx = [sample(xs,S) for _ in 1:T]
gs = gm.(sx)
ms = am.(sx)

AN = 50
anm=@animate for i in range(1, stop = AN)
    warpedx = minmaxreplace.(sx, i-1)
    gsp = gm.(warpedx)
    msp = am.(warpedx)

    p=boxplot([log.(gs),log.(gsp)], title="log-transformed", label=["X", "Y"])
    q=boxplot([log.(ms),log.(msp)], title="no-transform")
    Plots.plot(p, q, xticks=[],ylim=[4, 12], dpi=300, size=(800,400), ylabel="E[X]")
end
mkpath("figures")
g=gif(anm, joinpath("figures", "interactivelp.gif"), fps=5)
