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
module LogParadox
using Statistics
using SPECHT
using StatsBase
using Random

export gm, am, picki, smooth, minmaxreplace, to_entry, to_entries, generate_image, transform_steps, tf, reprand!, transform_steps, transform_steps_replace, tfsample, rep_min!, rep_max!, rep_minmax!, rep_rand!

"""
    gm(xs, base=exp(1))

    Computes the geometric mean for a vector of non-zero floating point values.
"""
function gm(xs, base = exp(1))
    exp(mean(log.(base, xs)))
end

"""
    am(xs)

    Arithmetic mean
"""
function am(xs)
    return mean(xs)
end

"""
    smooth(xs, step=2)

    Return a smoothed (windowed average with step size) version of xs
"""
function smooth(xs, step=2)
    N = length(xs)
    rs = zeros(N-step)
    for i in 1:N-step
        rs[i]=am(xs[i:i+step])
    end
    return rs
end

"""
    picki(N)

    Pick a random integer [0, N]
"""
function picki(N)
    return Int(round((rand()) * N))
end

"""
    warpx(xs, iterations=1)

    Replace the maximum and minimum of xs iteratively with the midpoint of [gm(xs), am(xs)]

    Returns a modified copy of xs
"""
function minmaxreplace(xs, iterations = 1)
    if iterations == 0
        return xs
    end
    _x = copy(xs)
    g = gm(_x)
    a = am(_x)
    G = (g + a) / 2
    for _ = 1:iterations
        m, i = findmin(_x)
        M, j = findmax(_x)
        _x[i] = G
        _x[j] = G
    end
    return _x
end


"""
    tf(xs, weight=1)

    Return the weighted midpoint of [GM, AM]
"""
function tf(xs, weight = 1)
    g = exp(mean(log.(xs)))
    m = mean(xs)
    p = (g + m) / 2
    return p * weight
end


"""
    transform_steps(xs, steps)

    Iteratively add the ID midpoint to xs
"""
function transform_steps(xs, steps)
    ys = copy(xs)
    mp = tf(ys)
    for _ = 1:steps
        push!(ys, mp)
        mp = tf(ys)
        # @info "New midpoint = mp"
    end
    return ys
end

"""
    transform_steps_replace(xs, steps, replacer, selector)

    For `steps` iterations, execute `replacer` on xs with values generated from `selector`
    Defaults to replace_minmax! and tfsample, a value sampled from the AM-GM interval.

    # Examples
    ```julia-repl
    ys = transform_steps_replace(randexp(1000).*100 .+ 1, 10)
    ```

    See also [`tfsample!`](@ref), [`rep_minmax!`](@ref), [`rep_rand!`](@ref).
"""
function transform_steps_replace(xs, steps, replacer = rep_minmax!, TF = tfsample)
    ys = copy(xs)
    mp = tf(ys)
    for _ = 1:steps
        mp = TF(ys)
        replacer(ys, mp)
    end
    return ys
end

"""
    tfsample(xs)

    Return a sampled value from the range [gm(xs), am(xs)]
"""
function tfsample(xs)
    g, m = gm(xs), am(xs)
    return g + rand() * (m - g)
end


"""
    rep_min!(xs, val)

    Replace the minimum value(s)
"""
function rep_min!(xs, val)
    xs[xs.==minimum(xs)] .= val
end

"""
    rep_max!(xs, val)

    Replace the maximum value(s)
"""
function rep_max!(xs, mp)
    xs[xs.==maximum(xs)] .= mp
end

"""
    rep_minmax!(xs, val)

    Replace the minimum and maximum value(s)
"""
function rep_minmax!(ys, mp)
    ys[ys.==maximum(ys)] .= mp
    ys[ys.==minimum(ys)] .= mp
end

"""
    rep_rand!(xs, val)

    Replace a randomly chosen value
"""
function rep_rand!(ys, mp)
    N = length(ys)
    ys[rand(1:N)] = mp
end

reprand! = rep_rand!



"""
    to_entries(samples, table)

    For each sample lookup in table what its frequency is.
"""
function to_entries(smpl, tbl)
    trs = [to_entry(s, tbl) for s in smpl]
    return trs
end


"""
    generate_image(counts, sizes, X, Y)

    Generate an image using SPECHT's in silico 2D generators

"""
function generate_image(counts, sizes, X, Y)
    GS=[]
    for a_c in keys(counts)
        σ = (sizes[a_c]^2)/5
        CT = counts[a_c]
        @info "Object size $a_c has freq $CT"
        cv = [[σ 0; 0 σ] for _ in 1:CT]
        rs = SPECHT.generate_rand_coordinates(X, Y, CT; offset=50)
        # GT = coordstogt([rs[i,:] for i in 1:N], X, Y)
        G = SPECHT.fastgaussian2d([rs[i,:] for i in 1:CT], cv, X, Y)
        push!(GS, G./maximum(G))
    end
    ima = GS[1] .+ GS[2] .+ GS[3] .+ GS[4]
end


"""
    to_entry(sample, table)

    Lookup entry in sorted frequency table
"""
function to_entry(s, tbl)
    @assert reverse(sort(tbl)) == tbl
    for (i, t) in enumerate(tbl)
        if s <= sum(tbl[1:i])
            return i
        end
    end
    return size(tbl)[1]
end

end
