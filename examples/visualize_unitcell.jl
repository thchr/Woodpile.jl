using Woodpile: Cylinder, Line, symmetrize
using Crystalline
using StaticArrays

sgnum = 222
kind = 4
sg = spacegroup(sgnum)
wps = wyckoffs(sgnum)

# ---------------------------------------------------------------------------------------- #
# some default choices for cylinder "seeds" for three kinds of space groups

if sgnum == 218
    # TODO: The center/axis choice is not good: produces equivalent symmetry to SG 223.
    rad = .25
    cntr = [0, 0, 0.5] # a (non-canonical) 6b Wyckoff position
    axis = [1, 0, 0] # [1,0,0] is not good; leads to sg 223 symmetry
elseif sgnum == 223
    rad = .2
    cntr = [0, 0, .5] # 2a Wyckoff position
    axis = [1, 0, 0]
elseif sgnum == 222
    # this transformation casts the SG to the same setting as SG 218 and SG 223
    # sg = transform.(sg, Ref(one(SMatrix{3,3,Float64})), Ref((@SVector [-1,-1,-1])/4))
    if kind == 1
        rad = .125
        cntr = wps[end-1]() # 6b Wyckoff position
        axis = [0, 1, 0] # flip to [1,0,0] for a "cubic-double-gyroid"-like structure
    elseif kind == 2
        rad = .2
        cntr = wps[end-1]()
        axis = [1, 0, 0] # "cubic-double-gyroid"-like structure
    elseif kind == 3
        rad = .1
        cntr = wps[end]()
        axis = [1, 1, 1] # "spoke-line"-like structure
    elseif kind == 4
        rad = .1
        cntr = wps[end-2]()
        axis = [1,0,0] # "cages"-like structure
    end
else
    error(DomainError(sgnum, "no default seed choices for `sgnum`; add your own"))
end
seed = Cylinder(Line(cntr, axis), rad)

# ---------------------------------------------------------------------------------------- #
# generate a symmetric woodpile structure based on `seed` and visualize it
using GLMakie

cs = symmetrize(sg, seed)
plot(cs, boxtol=1e-8)
