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
using Random

export gm, am, picki, minmaxreplace

"""
    gm(xs, base=exp(1))

    Computes the geometric mean for a vector of non-zero floating point values.
"""
function gm(xs, base=exp(1))
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
function minmaxreplace(xs, iterations=1)
    if iterations == 0
        return xs
    end
    _x = copy(xs)
    g=gm(_x)
    a=am(_x)
    G = (g + a)/2
    for _ in 1:iterations
        m, i = findmin(_x)
        M, j = findmax(_x)
        _x[i] = G
        _x[j] = G
    end
    return _x
end

end
