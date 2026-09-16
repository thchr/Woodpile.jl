# Inverse opal (SG 225, Fm-3m): close-packed air spheres on an FCC lattice, in a dielectric
# backbone. A useful reference structure - and a foil to the diamond one in
# `mpb-bandstructure-diamond.jl`, which differs only in the seed position and space group.

using PythonCall
using Brillouin, Crystalline, LinearAlgebra
using GLMakie
using Woodpile: Sphere, center, radius, symmetrize
using ProgressMeter

mp = pyimport("meep")
mpb = pyimport("meep.mpb")
mpb.verbosity = -1

## --------------------------------------------------------------------------------------- #
## sphere generation

sgnum = 225 # face-centered cubic Fm-3m
cntr = centering(sgnum, 3)

Rs′ = directbasis(sgnum; abclims = (1.0,1.0)) # conventional (cubic) basis, a = 1
Rs  = primitivize(Rs′, cntr)                  # primitive (FCC) basis
sg  = primitivize(spacegroup(sgnum))

# nearly close-packed: FCC spheres touch at a radius of 1/(2√2) (bond length 1/√2)
seed = Sphere([0,0,0], 1/(2sqrt(2)) * .96) # 4a Wyckoff position
sphs = symmetrize(sg, seed, Rs)
plot(sphs, Rs; inverted=true) # air spheres in a dielectric backbone

## --------------------------------------------------------------------------------------- #
# geometry

# k-vectors
kp  = irrfbz_path(sgnum, Rs′)
kvs = interpolate(kp, 55)

# meep geometry
# NB: MPB takes object coordinates in the (primitive) lattice basis, whereas Woodpile works
#     in Cartesian coordinates - hence the `latticize`. Radii are lengths in units of `a`,
#     so they need no conversion.
# NB: unlike PyCall, PythonCall does not auto-convert Julia vectors, so the `Vector3`s (and
#     the geometry list) must be built explicitly
m = mp.Medium(epsilon=1) # air spheres in an ε = 13 backbone (`default_material` below)
geometry = map(sphs) do sph
    mp.Sphere(center = mp.Vector3(latticize(center(sph), Rs)...),
              radius = radius(sph), material = m)
end
lattice = mp.Lattice(basis_size = mp.Vector3(norm.(Rs)...), # relative to conventional cell
                     basis1 = mp.Vector3(Rs[1]...),
                     basis2 = mp.Vector3(Rs[2]...),
                     basis3 = mp.Vector3(Rs[3]...))
ms = mpb.ModeSolver(
    num_bands        = 10,
    k_points         = [],
    geometry         = pylist(geometry),
    geometry_lattice = lattice,
    resolution       = 16,
    tolerance        = 1e-6,
    default_material = mp.Medium(epsilon=13),
)
ms.init_params(p = mp.NO_PARITY, reset_fields = true)
freqs = Matrix{Float64}(undef, length(kvs), pyconvert(Int, ms.num_bands))

@showprogress 0.1 for (i, kv) in enumerate(kvs)
    redirect_stdout(devnull) do
        ms.solve_kpoint(mp.Vector3(kv...))
    end
    freqs[i,:] = sort!(pyconvert(Vector{Float64}, ms.get_freqs()))
end

plot(kvs, freqs)
