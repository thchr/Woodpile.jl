using Woodpile: Cylinder, Ray, symmetrize
using GLMakie: plot
using Crystalline

sgnum = 223
sg = spacegroup(sgnum)
wps = get_wycks(sgnum)

# ---------------------------------------------------------------------------------------- #
# some default choices for cylinder "seeds" for three kinds of space groups

if sgnum ∈ (218, 223)
    rad = .25
    cntr = [0, 0, 0.5] # this corresponds to a 6b Wyckoff position
    axis = [1,0,0] # [0,1,0] is also good
    seed = Cylinder(Ray(cntr, axis), rad)
    #seed = [Cylinder(cntr, rad, axis), Cylinder(cntr, rad, [0,1,0])]
elseif sgnum == 222
    rad = .125
    cntr = wps[end-1]() # pick a 6b Wyckoff position
    axis = [0,1,0] # flip to [1,0,0] for a "cubic-double-gyroid"-like structure
    seed = Cylinder(Ray(cntr, axis), rad)
end

# ---------------------------------------------------------------------------------------- #
# generate a symmetric woodpile structure based on `seed` and visualize it

cs = symmetrize(sg, seed)
plot(cs)