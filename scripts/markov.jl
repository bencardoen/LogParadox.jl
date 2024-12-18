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
using Pkg; Pkg.activate(".")
using SPECHT
using Plots
using Random
using StatsBase
using Images, ImageView
using LogParadox
### This scripts generates in silico 2D images of fluorescent marked objects
### It uses SPECHT to generate the images, with the objects sampled from a Markov chain
Random.seed!(42)

X, Y = 512, 512

## Sizes of objects
sizes = [1, 3, 9, 27]
## Frequencies per celltype
freq_a = [300, 100, 30, 7]
freq_b = [240, 147, 30, 4]

Na = sum(freq_a)
Nb = sum(freq_b)

pa = freq_a ./ Na
pb = freq_b ./ Nb

S = 100
Random.seed!(42)

ac, bc, ima, imb, ga, gb = generate_images_from_markov_chains(pa, pb, sizes, sizes; X=X, Y=Y, S=S, matchstate=4)

@info ac[4] == bc[4]
# sizes
Images.save(joinpath("figures","a.tif"), N0f16.(imga./maximum(imga)))

Images.save(joinpath("figures","b.tif"), N0f16.(imgb./maximum(imgb)))
imshow(imgb)
imshow(imga)
