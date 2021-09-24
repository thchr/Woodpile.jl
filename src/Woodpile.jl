module Woodpile

using Crystalline: SymOperation, rotation, translation, isapproxin
using StaticArrays
using LinearAlgebra
using Requires: @require

import Base: minimum, maximum, *, ==, isequal, isapprox
# ---------------------------------------------------------------------------------------- #

struct Ray
    cntr::SVector{3,Float64}
    axis::SVector{3,Float64}
end
center(r::Ray) = r.cntr
axis(r::Ray) = r.axis

struct Cylinder
    ray::Ray
    radius::Float64
end
Cylinder(cntr::AbstractVector{<:Real}, axis::AbstractVector{<:Real}, radius::Real) = Cylinder(Ray(cntr, axis), radius)
ray(c::Cylinder) = c.ray
center(c::Cylinder) = center(ray(c))
axis(c::Cylinder)   = axis(ray(c))
radius(c::Cylinder) = c.radius

struct Box
    # axis aligned box
    cntr::SVector{3,Float64} # center
    width::SVector{3,Float64}
end
Box() = Box(SVector{3,Float64}(0.0, 0.0, 0.0))
Box(cntr::AbstractVector) = Box(cntr, SVector{3,Float64}(1.0, 1.0, 1.0))
center(b::Box) = b.cntr
width(b::Box)  = b.width
minimum(b::Box) = b.cntr .- width(b)./2
maximum(b::Box) = b.cntr .+ width(b)./2

struct Rect
    # a 2D rectangle with sides of length 1
    cntr::SVector{3,Float64}
    normal::SVector{3,Float64}
end

function facets(b::Box)
    normals = SVector{3, Float64}.(((1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)))
    cntrs = Ref(center(b)) .+ normals./2
    return Rect.(cntrs, normals)
end

# ---------------------------------------------------------------------------------------- #

function (*)(op::SymOperation{3}, c::Cylinder)
    cntr′ = rotation(op)*center(c) + translation(op)
    axis′ = rotation(op)*axis(c)
    axis′ = normalize(axis′)
    return Cylinder(Ray(cntr′, axis′), radius(c))
end
function (==)(c1::Cylinder, c2::Cylinder)
    # two cylinders' axes are the same if they are colinear
    iszero(axis(c1) × axis(c2)) || return false
    # two aligned cylinders' are the same if their centers lie on the line spanned
    # by their axis
    iszero((center(c1) - center(c2)) × axis(c1)) || return false    
    # two cylinders are the same if they have the same radius
    return radius(c1) == radius(c2)
end
function isapprox(c1::Cylinder, c2::Cylinder;
            atol::Real=1e-12, rtol::Real=atol>0 ? 0.0 : √eps())
    norm(axis(c1) × axis(c2)) ≤ atol || return false
    norm((center(c1) - center(c2)) × axis(c1)) ≤ atol || return false  
    return isapprox(radius(c1), radius(c2); atol=atol, rtol=rtol)
end

# ---------------------------------------------------------------------------------------- #

function distance_to_cylinder_axis(r::AbstractVector{<:Real}, c::Cylinder)
    # cf. https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
    tmp = r - center(c)
    return norm(tmp - dot(tmp, axis(c))*axis(c))
end
function signed_distance(r::AbstractVector{<:Real}, c::Cylinder)
    # evaluates the shortest distance from an input point `r` to a cylinder `c`: if `r` is
    # inside `c`, the returned value is positive; conversely, if `r` is outside `c`, the
    # returned value is negative
    return distance_to_cylinder_axis(r, c) - radius(c)
end
function signed_distance(r::AbstractVector{<:Real}, cs::AbstractVector{Cylinder})
    # evaluates the minimum signed distance from `r` to any number of cylinders `cs`
    d = Inf
    for c in cs
        d = min(d, signed_distance(r, c))
    end
    return d
end

function intersects(b::Box, r::Ray)
    # Smits' algorithm from https://dl-acm-org/doi/pdf/10.1145/1198555.1198748
    b_min, b_max = minimum(b), maximum(b)
    r_cntr = center(r)
    r_axis = axis(r)
    if r_axis[1] ≥ 0
        tmin = (b_min[1] - r_cntr[1]) / r_axis[1]
        tmax = (b_max[1] - r_cntr[1]) / r_axis[1]
    else
        tmin = (b_max[1] - r_cntr[1]) / r_axis[1]
        tmax = (b_min[1] - r_cntr[1]) / r_axis[1]
    end

    if r_axis[2] ≥ 0
        tymin = (b_min[2] - r_cntr[2]) / r_axis[2]
        tymax = (b_max[2] - r_cntr[2]) / r_axis[2]
    else
        tymin = (b_max[2] - r_cntr[2]) / r_axis[2]
        tymax = (b_min[2] - r_cntr[2]) / r_axis[2]
    end

    (tmin > tymax || tymin > tmax) && return (bool=false, tmin=NaN, tmax=NaN)
    tmin = max(tmin, tymin)
    tmax = min(tmax, tymax)
    
    if r_axis[3] ≥ 0
        tzmin = (b_min[3] - r_cntr[3]) / r_axis[3]
        tzmax = (b_max[3] - r_cntr[3]) / r_axis[3]
    else
        tzmin = (b_max[3] - r_cntr[3]) / r_axis[3]
        tzmax = (b_min[3] - r_cntr[3]) / r_axis[3]
    end

    (tmin > tzmax || tzmin > tmax) && return (bool=false, tmin=NaN, tmax=NaN)
    tmin = max(tzmin, tmin)
    tmax = min(tmax, tzmax)
    return (bool=true, tmin=tmin, tmax=tmax)
end
function might_intersects(b::Box, c::Cylinder)
    # poor man's intersection check: if a cylinder's center doesn't penetrate the faces
    # of a box which is padded by the cylinder's radius, then it cannot possibly intersect
    # the box itself (it might still not intersect the box, however).
    # to do better than this, we should check if the cylinder itself ever intersects with
    # the edge-segments of the box (but that is tedious... so we just make do with this)
    b′ = Box(center(b), width(b) .+ radius(c)*2)
    return intersects(b′, ray(c))
end
might_intersects(c::Cylinder, b::Box=Box()) = might_intersects(b, c)

# ---------------------------------------------------------------------------------------- #

const NEIGHBOR_TRANSLATIONS = 
    [SVector(I[1],I[2],I[3]) for I in vec(CartesianIndices((-1:1, -1:1, -1:1))) if !iszero(I)]

function symmetrize(ops::AbstractVector{SymOperation{3}}, cs::AbstractVector{Cylinder};
            add_neighbors::Bool=true)
    if add_neighbors
        ops = vcat(ops, SymOperation{3}.(NEIGHBOR_TRANSLATIONS))
    end
    cs′ = deepcopy(cs)
    for op in ops
        for c in cs
            c′ = op * c
            if !isapproxin(c′, cs′) && might_intersects(c′).bool
                push!(cs′, c′)
            end
        end
    end

    if length(cs′) == length(cs)
        # no extra cylinders were added, we converged       
        # however, before returning we do the following check: in case the initial seed
        # point was outside the unit cell, there might be an element in `cs′` that really
        # doesn't intersect the unitcell; no point in keeping it around - filter it out here
    
        return filter!(c->might_intersects(c).bool, cs′) # no extra cylinders were added, we converged
    else
        return symmetrize(ops, cs′; add_neighbors=false)
    end
end
symmetrize(ops::AbstractVector{SymOperation{3}}, c::Cylinder; kwargs...) = symmetrize(ops, [c]; kwargs...)

# ---------------------------------------------------------------------------------------- #

function __init__()
    @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" begin
        @eval using GeometryBasics: Rect, Vec, Mesh
        @eval using Meshing
        @eval import GLMakie: plot
        function plot(cs::AbstractVector{Cylinder}; style::Symbol=:merged)
            style = :merged
            
            f = Figure()
            ax = Axis3(f; aspect=:data, limits=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5))
            
            plotstyle = :merged # :merged or :individual
            if style == :individual
                surfs = Mesh.(Base.Fix2.(Ref(signed_distance), cs),
                            Ref(Rect(Vec(-.5,-.5,-.5), Vec(1,1,1))),
                            Ref(MarchingCubes(iso=0)), samples=(50,50,50))
                foreach(surfs) do surf
                    GLMakie.mesh!(surf, colormap=:hot, shading=true)
                end
            
            elseif style == :merged
                surf = Mesh(Base.Fix2(signed_distance, cs),
                            Rect(Vec(-.5,-.5,-.5), Vec(1,1,1)),
                            MarchingCubes(iso=0), samples=(50,50,50))
                GLMakie.mesh!(surf, colormap=:hot, shading=true)
            end
            xlims!(ax, (-0.5,0.5)); ylims!(ax, (-0.5,0.5)); zlims!(ax, (-0.5,0.5))

            return f
        end
    end
end
# ---------------------------------------------------------------------------------------- #
end # Module