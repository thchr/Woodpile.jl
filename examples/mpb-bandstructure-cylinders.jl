using PyCall
using Brillouin, Crystalline, LinearAlgebra
using Woodpile: Cylinder, Line, symmetrize

mp = pyimport("meep")
mpb = pyimport("meep.mpb")
mpb.verbosity = -1

## --------------------------------------------------------------------------------------- #
## cylinder generation
function create_cylinder_unitcell(sgnum, Rs, rad=nothing, kind=1)
    sg = spacegroup(sgnum)
    wps = wyckoffs(sgnum)
    if sgnum == 218
        # TODO: The center/axis choice is not good: produces equivalent symmetry to SG 223.
        rad = something(rad, .25)
        cntr = [0, 0, 0.5] # a (non-canonical) 6b Wyckoff position
        axis = [1, 0, 0] # [1,0,0] is not good; leads to sg 223 symmetry
    elseif sgnum == 223
        rad = something(rad, .2)
        cntr = [0, 0, .5] # 2a Wyckoff position
        axis = [1, 0, 0]
    elseif sgnum == 222
        # this transformation casts the SG to the same setting as SG 218 and SG 223
        #     transform.(sg, Ref(one(SMatrix{3,3,Float64})), Ref((@SVector [-1,-1,-1])/4))
        if kind == 1
            rad = something(rad, .2)
            cntr = wps[end-1]() # 6b Wyckoff position
            axis = [0, 1, 0] # "doubly-connected tri-crosses" structure
        elseif kind == 2
            rad = something(rad, .2)
            cntr = wps[end-1]()
            axis = [1, 0, 0] # "cubic-double-gyroid" structure
        elseif kind == 3
            rad = something(rad, .2)
            cntr = wps[end]()
            axis = [1, 1, 1] # "spoke-line" structure
        elseif kind == 4
            rad = something(rad, .2)
            cntr = wps[end-2]()
            axis = [1,0,0] # "cages" structure
        end
    else
        error(DomainError(sgnum, "no default seed choices for `sgnum`; add your own"))
    end

    seed = Cylinder(Line(cntr, axis), rad)
    return symmetrize(sg, seed, Rs)
end

## --------------------------------------------------------------------------------------- #
# geometry

sgnum = 222
Rs′ = directbasis(sgnum)
Rs  = primitivize(Rs′, centering(sgnum))

cyls = create_cylinder_unitcell(sgnum, Rs, nothing, 3)

# k-vectors
kp  = irrfbz_path(sgnum, Rs′)
kvs = interpolate(kp, 35)

# meep geometry
m = mp.Medium(epsilon=13)
geometry = map(cyls) do cyl
    r = cyl.radius
    c = cyl.line.cntr
    a = cyl.line.axis
    mp.Cylinder(center=c, radius=r, material=m, height=mp.inf, axis=mp.Vector3(a...))
end

lattice = mp.Lattice(basis_size=norm.(Rs), # take units relative to conventional unit cell
                     basis1 = Rs[1], basis2 = Rs[2], basis3 = Rs[3])
ms = mpb.ModeSolver(
    num_bands        = 10,
    k_points         = [],
    geometry         = geometry,
    geometry_lattice = lattice,
    resolution       = 16,
    tolerance        = 1e-6,
    )

ms.init_params(p = mp.NO_PARITY, reset_fields = true)
freqs = Matrix{Float64}(undef, ms.num_bands, length(kvs))

for (i, kv) in enumerate(kvs)
    ms.solve_kpoint(mp.Vector3(kv...))
    freqs[:,i] = ms.get_freqs()
end

## --------------------------------------------------------------------------------------- #
# UnicodePlots
# kdists = Brillouin.cumdists(kvs)
# p = lineplot(kdists, freqs[1,:], color=2, xlim=extrema(kdists), ylim=extrema(freqs), 
#              width=displaysize(stdout)[2]-20, height=displaysize(stdout)[1]-6)
# for b in 2:ms.num_bands
#     lineplot!(p, kdists, freqs[b,:], color=2)
# end
# display(p)

## --------------------------------------------------------------------------------------- #
import PlotlyJS
p = PlotlyJS.plot(kvs, freqs')
display(p)

## --------------------------------------------------------------------------------------- #

using MPBUtils
lgirsd = lgirreps(sgnum)
lgirs = realify(lgirsd["R"])
lg = group(first(lgirs))
symeigs = [Vector{ComplexF64}(undef, length(lg)) for _ in 1:ms.num_bands] # [band][opidx]
kv = mp.Vector3(position(lg)()...)

#ms.num_bands  = 6
#ms.resolution = 32
ms.init_params(p = mp.NO_PARITY, reset_fields = true)
ms.solve_kpoint(kv)

for (opidx, op) in enumerate(lg)
    W = mp.Matrix(op.rotation[:,1], op.rotation[:,2], op.rotation[:,3])
    w = mp.Vector3(op.translation...)
    symeigs′ = ms.compute_symmetries(W, w)
    for band in 1:ms.num_bands
        symeigs[band][opidx] = symeigs′[band]
    end
end
MPBUtils.find_individual_multiplicities(Dict(""=>symeigs), Dict(""=>lgirs))[""]

## --------------------------------------------------------------------------------------- #
# Berry phases
function spherical_loops(cntr, rad, Nϕ, Nθ)
    ϕs = range(0, 2π, Nϕ)
    θs = range(π, 0, Nθ+1)
    θs = θs[1:end-1] .- π/(2Nθ) # don't want the poles (zero-measure loops) ...
    return [[[cntr[1] + rad*cos(ϕ)*sin(θ),
              cntr[2] + rad*sin(ϕ)*sin(θ),
              cntr[3] + rad*cos(θ)        ] for ϕ in ϕs] for θ in θs]
end

Nϕ, Nθ = 10, 30
k_lists = spherical_loops([1,1,1]/2, .1, Nϕ, Nθ)
k_lists′ = [[mp.Vector3(k[1], k[2], k[3]) for k in k_list] for k_list in k_lists];

@pyinclude("/home/tchr/local/playground/mpb_berry_phase.py")

py"berry_phase($ms, $(k_lists′[1]), [[1,2,3]])"

phases = map(k_lists′) do k_list′
    py"berry_phase($ms, $k_list′, [[1,2,3], [4,5,6], [1,2,3,4,5,6]])"
end

using UnicodePlots
θs = range(π, 0, Nθ+1)[1:end-1] .- π/(2*(Nθ-1))
lineplot(θs, sum.(getindex.(phases, 1)), ylim=[-π, π]/1.1, xlim=[0,π])
