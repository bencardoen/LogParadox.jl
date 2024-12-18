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
# Copyright 2022-4, Ben Cardoen
module LogParadox
using Statistics
using SPECHT
using StatsBase
using Combinatorics
using Random
using DataFrames

export gm, am, picki, smooth, ID, check_paradox, generate_images_from_markov_chains, check_dataframe, minmaxreplace, to_entry, to_entries, generate_image, transform_steps, tf, reprand!, transform_steps, transform_steps_replace, tfsample, rep_min!, rep_max!, rep_minmax!, rep_rand!

"""
    gm(xs, base=exp(1))

    Computes the geometric mean for a vector of non-zero floating point values.
"""
function gm(xs, base = exp(1))
    exp(mean(log.(base, xs)))
end

"""
    check_paradox(xs, ys)
    Return the means (geometric and arithmetic) for x, y and if there is a paradoxical ordering.
"""
function check_paradox(xs::AbstractVector{T}, ys::AbstractVector{T}) where {T<:Number}
    gx, ax = gm(xs), am(xs)
    gy, ay = gm(ys), am(ys)
    p = false
    if gx < gy
        p = gx < gy && ay < ax
    else
        p = gy < gx && ax < ay
    end
    return gx, gy, ax, ay, p
end 

"""
    check_dataframe(df, labelcolumn="label", skipcols=nothing)
    
    Check if a dataframe contains a paradoxical comparison. 
    The dataframe should have N rows where each column, except `labelcolumn` is a feature
"""
function check_dataframe(df, labelcolumn="label", skipcols=nothing)
    columns = names(df)
    if ! (labelcolumn in columns)
        @error "No label"
        throw(ArgumentError("No label to check"))
    end
    labels = unique(df[:, labelcolumn])
    features = [c for c in columns if c != labelcolumn]
    @info "Have $labels different labels for $(size(df, 1)) rows and $(length(features)) features"
    label_combinations = combinations(labels, 2)
    # @info label_combinations |> collect
    results = DataFrame([String[],String[],String[], Float64[], Float64[], Float64[], Float64[], Bool[]], ["feature", "label_1", "label_2", "gx", "gy", "ax", "ay", "paradox?"])
    for f in features
        if !isnothing(skipcols) && f in skipcols
                continue
        end
        # @info f
        for (lx, ly) in label_combinations
            # @info lx, ly
            dfx = df[df[:,labelcolumn] .== lx, f]
            dfy = df[df[:,labelcolumn] .== ly, f]
            # @info typeof(dfx)
            gx, gy, ax, ay, p = check_paradox(dfx, dfy)
            push!(results, [f, "$lx", "$ly", gx, gy, ax, ay, p])
        end
    end        
    return results
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
    ID(xs)

    Intermean distance (am(xs)-gm(xs))
"""
function ID(xs)
    am(xs) - gm(xs)
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
    Counts is Dict of 1:K objects, corresponding with indices in `sizes`.
    A size of `2` will correspond to a Gaussian PSF at a random location, K times, of size (2^2)/denom
    X, Y are image dimensions.

"""
function generate_image(counts, sizes, X, Y; gt=false, denom=5, offset=50)
    GS=[]
    coords = Dict()
    for a_c in keys(counts)
        σ = (sizes[a_c]^2)/denom
        CT = counts[a_c]
        @info "Object size $a_c has freq $CT"
        cv = [[σ 0; 0 σ] for _ in 1:CT]
        rs = SPECHT.generate_rand_coordinates(X, Y, CT; offset=offset)
        GT = coordstogt([rs[i,:] for i in 1:CT], X, Y)
        G = SPECHT.fastgaussian2d([rs[i,:] for i in 1:CT], cv, X, Y)
        push!(GS, G./maximum(G))
        coords[a_c] = GT
    end
    ima = GS[1] .+ GS[2] .+ GS[3] .+ GS[4]
    if gt
        return ima, coords
    else
        return ima
    end
end

"""
    Create images from a Markov chain distribution with two different chains and sizes.
    Matchstates is a control to generate until a frequency of objects across a, b, is the same.
    pa, pb should be probability vectors of length N
    sizes_{a,b} are the object sizes, aligned with pa, pb.
    pb[i] is the probability of observing sizes of sizes_b[i].
    X, Y are image dimensions (>128)
    S is the number of total objects to generate.
    Return the frequencies and images.
"""
function generate_images_from_markov_chains(pa, pb, sizes_a, sizes_b; X=512, Y=512, S=100, matchstate=nothing)
    a_counts = nothing
    b_counts = nothing
    while true
        rs = rand(S)
        entries_a = to_entries(rs, pa)
        entries_b = to_entries(rs, pb)
        a_counts = countmap(entries_a)
        b_counts = countmap(entries_b)
        if isnothing(matchstate)
            break
        else
            if length(a_counts)==length(b_counts)
                if a_counts[matchstate] == b_counts[matchstate]
                    break
                end
            end
        end
    end
    imga, ga = generate_image(a_counts, sizes_a, X, Y; gt=true)
    imgb, gb = generate_image(b_counts, sizes_b, X, Y; gt=true)
    return a_counts, b_counts, imga, imgb, ga, gb
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
