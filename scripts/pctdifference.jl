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
_f = transform_steps
F = reprand!
_f = (x, y) -> transform_steps_replace(x, y, F)
S = Int.(length(es)/2)
Random.seed!(42)
steps = [mean(_f(es, i)) for i in 1:S]
ds = [steps[1]/steps[i] for i in 1:S]
gsteps = [gm(_f(es, i)) for i in 1:S]
dsg = [gsteps[1]/gsteps[i] for i in 1:S]
z = Plots.plot(100 .* (1 .- ds), alpha=.1, label="δ μ+")
Plots.plot!(100 .*(1 .- dsg), alpha=.1, label="δ μ*", xlabel="Number of elements replaced", ylabel="Differential of mean (%)")
Plots.plot!(smooth(100 .*(1 .- dsg),10), alpha=1, label=false, color=:orange, xlabel="Number of elements replaced", ylabel="Differential of mean (%)")
Plots.plot!(smooth(100 .*(1 .- ds),10), alpha=1, label=false, color=:blue, xlabel="Number of elements replaced", ylabel="Differential of mean (%)")

# Plots.plot!(100 .* (abs.(1 .- dsg./ds)), label="Ratio: maximizing difference", xlabel="Number of elements altered", ylabel="% change in mean")
Plots.plot(z, dpi=300, size=(800, 600))
Plots.savefig(joinpath("figures", "differentialrand.pdf"))
# Plots.plot(steps)
