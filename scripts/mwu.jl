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
using Random, Distributions, Plots, HypothesisTests
using LogParadox
Random.seed!(42)
es=randexp(1000)*1000 .+ 10
NPX=200
pxs = zeros(NPX)
pys = zeros(NPX)
SNX = 50
SAM = 200
for i in 1:NPX
    ES=transform_steps_replace(es, i)
    gps = [gm(sample(ES, SAM)) for _ in 1:SNX]
    mps = [mean(sample(ES, SAM)) for _ in 1:SNX]
    gs = [gm(sample(es, SAM)) for _ in 1:SNX]
    ms = [mean(sample(es, SAM)) for _ in 1:SNX]
    X = HypothesisTests.MannWhitneyUTest(ms, mps)
    pxs[i] = pvalue(X)
    Y = HypothesisTests.MannWhitneyUTest(gs, gps)
    pys[i] = pvalue(Y)
end

Xs = [1:SNX]
z=Plots.plot(smooth(log.(pxs), 10),label="Arithmetic mean difference")
Plots.plot!(smooth(log.(pys),10), label="Geometric mean difference")
Plots.scatter!(log.(pxs),  alpha=.1, label=false, color="blue")
Plots.scatter!(log.(pys), alpha=.1, label=false, color="red")
# Plots.scatter!(smooth(log.(pys),10), label=false)
Plots.plot!(log.([0.05 for _ in 1:NPX]), label="*")
Plots.plot!(log.([0.01 for _ in 1:NPX]), label="**")
Plots.plot!(log.([0.001 for _ in 1:NPX]), label="***")
Plots.plot(z, ylim=(-40, 0), dpi=300, size=(800, 600), title="Divergence under resampled (N=$SNX, S=$SAM) means",  xlabel="Elements replaced (x/2000)", ylabel="P-value Mann Whitney U test (log)")
# z=Plots.plot(sort(log.(es)), sort(es), markersize=1, xlabel="Log(A)", ylabel="A", label="Sorted A")
Plots.savefig("pvals$(SNX)_$(SAM).pdf")
