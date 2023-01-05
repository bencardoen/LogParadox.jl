using Images, Distributions, Plots, LogParadox, GR


m = 1
xs = 1:25 |> collect

# f(m, M) = am([m, M]) - gm([m, M])

ys = [ID([m, M]) for M in xs]

plot!(xs, ys, label="d(1, M)")
p=Plots.plot(xs, xs, label="M=maximum(X)")
Plots.plot!(xs, sqrt.(xs), label="√(1*M)")
Plots.plot!(xs, (1 .+ xs)./2, label="(1+M)/2")
Plots.plot(p, dpi=600, size=(800, 600), legend=:topleft)


gr()
x=range(1,stop=500,length=100)
y=range(1,stop=500,length=100)
f(x,y) = log(10, 1+ID([x,y]))

p=Plots.plot(x,y,f,st=:surface,camera=(-30, 35), alpha=.5,  xlabel="Minimum X", ylabel="Maximum X", zlabel="d(m,M)",c=cgrad([:green,:red]),size=(600, 600), dpi=100)
Plots.savefig(p, joinpath("figures","mesh.svg"))
