using Woodpile: Cylinder, Ray, symmetrize
using Crystalline

sgnum = 223
kind = 1
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
    if kind == 1
        rad = .125
        cntr = wps[end-1]() # 6b Wyckoff position
        axis = [0, 1, 0] # flip to [1,0,0] for a "cubic-double-gyroid"-like structure
    elseif kind == 2
        # this transformation casts the SG to the same setting as SG 218 and SG 223
        sg = transform.(sg, Ref(SMatrix{3,3,Float64}(I(3))), Ref((@SVector [-1,-1,-1])/4))
        rad = .2
        cntr = [0, 0, 0.5]
        axis = [1, 0, 0]
    end
else
    error(DomainError(sgnum, "no default seed choices for `sgnum`; add your own"))
end
seed = Cylinder(Ray(cntr, axis), rad)

# ---------------------------------------------------------------------------------------- #
# generate a symmetric woodpile structure based on `seed` and visualize it
using GLMakie: plot

cs = symmetrize(sg, seed)
plot(cs)