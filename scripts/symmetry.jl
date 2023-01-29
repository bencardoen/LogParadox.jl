using Plots, GR, Statistics
using Distributions
using StatsPlots
using LogParadox

idx = []
ax = []
gx = []
for i in 100:1000
    d=Normal(i, 10)
    r=rand(d, 100)
    push!(idx, am(r)-gm(r))
    push!(ax, am(r))
    push!(gx, gm(r))
end
xa = []
ga = []
for i in 10:100
    X = rand(Normal(i, 2), 100)
    X[1] = sqrt(i)
    X[end] = i^2
    push!(xa, am(X))
    push!(ga, gm(X))
end

Plots.plot(idx)

p=Plots.scatter(xa .- ga)
Plots.plot(p, dpi=90, xlabel="sample(N(μ = x, σ = 2), 100) || [√(x), x^2] ", ylabel="ID(X)", title="Symmetric distribution's sensitive to paradox")
Plots.savefig(joinpath("figures", "symmetric.png"))
