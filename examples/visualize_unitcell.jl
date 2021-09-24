
using Woodpile: Cylinder, Ray, symmetrize
using GLMakie: plot


## 
sgnum = 223
sg = spacegroup(sgnum)
wps = get_wycks(sgnum)

##
if sgnum ∈ (218, 223)
    rad = .25
    cntr = [0, 0, 0.5] # this corresponds to a 6b Wyckoff position
    axis = [1,0,0] # [0,1,0] is also good
    seed = Cylinder(Ray(cntr, axis), rad)
    #seed = [Cylinder(cntr, rad, axis), Cylinder(cntr, rad, [0,1,0])]
elseif sgnum == 222
    rad = .125
    cntr = wps[end-1]() # pick a 6b Wyckoff position
    axis = [0,1,0]
    seed = Cylinder(Ray(cntr, axis), rad)
end

##
cs = Woodpile.symmetrize(sg, seed; add_neighbors=true)
plot(cs)