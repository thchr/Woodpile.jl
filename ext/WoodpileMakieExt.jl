module WoodpileMakieExt

# ---------------------------------------------------------------------------------------- #

import Makie: plot, plot!
using Makie: Makie, Figure, Axis3, mesh!, lines!, FigureAxisPlot, GeometryBasics,
             xlims!, ylims!, zlims!
using Meshing: isosurface
using StaticArrays: SVector

using Brillouin: Brillouin, Cell, basis, in_wignerseitz, setting, vertices
using Bravais: DirectBasis, cartesianize

using Woodpile: signed_distance, Cylinder, Sphere, facets

# ---------------------------------------------------------------------------------------- #

function plot(
    cs::AbstractVector{<:Union{Cylinder, Sphere}},
    boundary::Union{AbstractVector{<:AbstractVector{<:Number}}, DirectBasis{3}, Cell{3}, Nothing} = nothing
    ; # keyword arguments
    style::Symbol=:merged,
    samples::Union{NTuple{3, Int}, Int}=50,
    axis::NamedTuple=NamedTuple(),
    figure::NamedTuple=NamedTuple(),
    plot_kws...,
)
    f = Figure(; figure...)
    ax = _default_bare_axis!(f, Val(3); axis=axis)
    p = plot!(ax, cs, boundary; style=style, samples=samples, plot_kws...)

    return FigureAxisPlot(f, ax, p)
end

# ---------------------------------------------------------------------------------------- #
# Trapezoidal boundary

function plot!(
    ax::Axis3,
    cs::AbstractVector{<:Union{Cylinder, Sphere}},
    boundary::Union{AbstractVector{<:AbstractVector{<:Number}}, DirectBasis{3}, Nothing} = nothing
    ; # keyword arguments
    style::Symbol=:merged,
    samples::Union{NTuple{3, Int}, Int}=50,
    plot_kws...,
)
    Rs = if boundary isa DirectBasis{3}
        boundary
    elseif boundary isa AbstractVector{<:AbstractVector{<:Number}}
        DirectBasis{3, Float64}(boundary)
    else
        DirectBasis{3, Float64}(SVector(1,0,0), SVector(0,1,0), SVector(0,0,1))
    end :: DirectBasis{3}
    
    # 3D unit cell
    if boundary !== nothing
        fs = facets(Rs)
        segments = Vector{SVector{3,Float64}}()
        for (i, f) in enumerate(fs)
            append!(segments, f)
            push!(segments, f[1])
            i == length(fs) || push!(segments, SVector(NaN, NaN, NaN))
        end
        lines!(ax, segments; color=:black, linewidth=0.5)

        xlims!(ax, minimum(f->minimum(r->r[1], f), fs), maximum(f->maximum(r->r[1], f), fs))
        ylims!(ax, minimum(f->minimum(r->r[2], f), fs), maximum(f->maximum(r->r[2], f), fs))
        zlims!(ax, minimum(f->minimum(r->r[3], f), fs), maximum(f->maximum(r->r[3], f), fs))
    else
        xlims!(ax, -0.5, 0.5); ylims!(ax, -0.5, 0.5); zlims!(ax, -0.5, 0.5)
    end

    # show cylinders / spheres, cropped to 3D unit cell bounding box
    plot_opts = merge((; color=:gray, transparency=false), plot_kws)

    samples isa Int && (samples = (samples, samples, samples))
    δ = 1e-8
    xs = range(-0.5-δ, 0.5+δ, length=samples[1])
    ys = range(-0.5-δ, 0.5+δ, length=samples[2])
    zs = range(-0.5-δ, 0.5+δ, length=samples[3])
    #= TODO: use below if https://github.com/JuliaGeometry/Meshing.jl/issues/104 is fixed
    δ = 1e-8
    xs = vcat(-0.5-δ, range(-0.5+δ, 0.5-δ, length=samples[1]), 0.5+δ)
    ys = vcat(-0.5-δ, range(-0.5+δ, 0.5-δ, length=samples[2]), 0.5+δ)
    zs = vcat(-0.5-δ, range(-0.5+δ, 0.5-δ, length=samples[3]), 0.5+δ)
    =#
    p = if style == :individual
        local _p
        for c in cs
            V = [signed_distance_if_inside(x, y, z, c, Rs) for x in xs, y in ys, z in zs]
            vs, fs_idxs = isosurface(V, xs, ys, zs)
            fs_idxs = map(i -> GeometryBasics.TriangleFace(i[1], i[2], i[3]), fs_idxs)
            _p = mesh!(cartesianize.(SVector.(vs), Ref(Rs)), fs_idxs; plot_opts...)
        end
        _p
    elseif style == :merged
        V = [signed_distance_if_inside(x, y, z, cs, Rs) for x in xs, y in ys, z in zs]
        vs, _fs_idxs = isosurface(V, xs, ys, zs)
        fs_idxs = map(i -> GeometryBasics.TriangleFace(i[1], i[2], i[3]), _fs_idxs)
        mesh!(cartesianize.(SVector.(vs), Ref(Rs)), fs_idxs; plot_opts...)
    else
        error("unsupported style $style")
    end

    return p
end

function inside_box(x, y, z)
    (x < -0.5 || x > 0.5) && return false
    (y < -0.5 || y > 0.5) && return false
    (z < -0.5 || z > 0.5) && return false
    return true
end

function signed_distance_if_inside(
    x::Real, y::Real, z::Real, 
    c::Union{Cylinder, AbstractVector{Cylinder}},
    Rs::DirectBasis{3}
)
    # NB: we assume `x`, `y`, and `z` to be given in _lattice_ coordinates here, unlike
    #     in the `uc::Cell{3}` method variant below
    inside_box(x, y, z) || return 1e20 # outside box
    rc = cartesianize(SVector(x,y,z), Rs)
    return signed_distance(rc, c) # inside box: return distance to cylinder(s)
end

# ---------------------------------------------------------------------------------------- #
# Wigner-Seitz cell boundary

function plot!(
    ax::Axis3,
    cs::AbstractVector{<:Union{Cylinder, Sphere}},
    uc::Cell{3}
    ; # keyword arguments
    style::Symbol=:merged,
    samples::Union{NTuple{3, Int}, Int}=50,
    plot_kws...,
)
    if setting(uc) == Brillouin.LATTICE
        uc = cartesianize(uc) # always only work with cartesian unit cells
    end
    
    # 3D unit cell
    xmin, xmax, ymin, ymax, zmin, zmax = unitcell_bounding_box(uc)
    plot!(ax, uc)
    xlims!(ax, xmin, xmax); ylims!(ax, ymin, ymax); zlims!(ax, zmin, zmax)

    # show cylinders / spheres, cropped to 3D unit cell bounding box
    plot_opts = merge((; color=:gray, transparency=false), plot_kws)

    samples isa Int && (samples = (samples, samples, samples))
    δ = 1e-8
    xs = range(xmin-δ, xmax+δ, length=samples[1]) # NB: Cartesian coordinates, not latticized
    ys = range(ymin-δ, ymax+δ, length=samples[2])
    zs = range(zmin-δ, zmax+δ, length=samples[3])
    #= TODO: use if https://github.com/JuliaGeometry/Meshing.jl/issues/104 is fixed
    xs = vcat(xmin-δ, range(xmin+δ, xmax-δ, length=samples[1]), xmax+δ)
    ys = vcat(ymin-δ, range(ymin+δ, zmax-δ, length=samples[2]), zmax+δ)
    zs = vcat(zmin-δ, range(zmin+δ, zmax-δ, length=samples[3]), zmax+δ)
    =#
    # NB: Unlike the `Rs` `plot!(...)`-variant above, here, the coordinates are cartesian
    #     from the get-go; i.e., we don't need to go back & forth between coordinate systems
    p = if style == :individual
        local _p
        for c in cs
            V = [signed_distance_if_inside(x, y, z, c, uc) for x in xs, y in ys, z in zs]
            vs, fs_idxs = isosurface(V, xs, ys, zs)
            fs_idxs = map(i -> GeometryBasics.TriangleFace(i[1], i[2], i[3]), fs_idxs)
            _p = mesh!(vs, fs_idxs; plot_opts...)
        end
        _p
    elseif style == :merged
        V = [signed_distance_if_inside(x, y, z, cs, uc) for x in xs, y in ys, z in zs]
        vs, _fs_idxs = isosurface(V, xs, ys, zs)
        fs_idxs = map(i -> GeometryBasics.TriangleFace(i[1], i[2], i[3]), _fs_idxs)
        mesh!(vs, fs_idxs; plot_opts...)
    else
        error("unsupported style $style")
    end

    return p
end

function signed_distance_if_inside(
    x::Real, y::Real, z::Real, 
    c::Union{Cylinder, AbstractVector{Cylinder}},
    uc::Cell{3}
)
    # NB: we assume `x`, `y`, and `z` to be given in _cartesian_ coordinates here, unlike
    #     in the `Rs::DirectBasis{3}` method-variant above
    rᶜ = SVector(x, y, z)
    in_wignerseitz(rᶜ, uc) || return 1e20 # outside wigner-seitz cell
    return signed_distance(rᶜ, c) # inside wigner-seitz cell: return distance to cylinder(s)
end

function unitcell_bounding_box(uc::Cell{3})
    if setting(uc) == Brillouin.LATTICE
        error("input must be in CARTESIAN setting: pass through `cartesianize!(uc)` first")
    end
    vs = vertices(uc)
    xmin, xmax = extrema(v->v[1], vs)
    ymin, ymax = extrema(v->v[2], vs)
    zmin, zmax = extrema(v->v[3], vs)

    return (xmin, xmax, ymin, ymax, zmin, zmax)
end

# ---------------------------------------------------------------------------------------- #

# copied from Brillouin/BrillouinMakieExt/wignerseitz.jl
function _default_bare_axis!(f, ::Val{3}; hideaxis::Bool=true, axis=NamedTuple())
    f[1,1] = ax = Axis3(f;
        aspect=:data,
        viewmode=:free,
        axis...)
    if hideaxis
        Makie.hidedecorations!(ax); ax.protrusions[] = 0 # cf. https://github.com/MakieOrg/Makie.jl/issues/2259
        Makie.hidespines!(ax)
    end
    return ax
end

end # module WoodpileMakieExt