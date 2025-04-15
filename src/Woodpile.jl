module Woodpile

# ---------------------------------------------------------------------------------------- #

using Crystalline: SymOperation, rotation, translation, isapproxin, @S_str, compose
using Bravais: DirectBasis, cartesianize, transform, latticize
using Brillouin: Cell, Brillouin, setting

using StaticArrays: SVector, StaticVector
using LinearAlgebra: norm, normalize, dot, ×

import Base: minimum, maximum, *, ==, isapprox
import Bravais: cartesianize, latticize

# ---------------------------------------------------------------------------------------- #

export Cylinder, Line, Sphere
export center, axis, radius
export symmetrize


## --------------------------------------------------------------------------------------- #

struct Line
    cntr::SVector{3,Float64}
    axis::SVector{3,Float64}
    function Line(cntr::SVector{3,Float64}, axis::SVector{3,Float64})
        l = norm(axis)
        iszero(l) && error("Line axis cannot be a zero vector")
        new(cntr, isone(l) ? axis : axis./l)
    end
end
Line(cntr, axis) = Line(convert(SVector{3,Float64}, cntr), convert(SVector{3,Float64}, axis))
center(l::Line)  = l.cntr
axis(l::Line)    = l.axis
cartesianize(l::Line, basis::DirectBasis{3}) = Line(center(l)'*basis, axis(l)'*basis)
latticize(l::Line, basis::DirectBasis{3})    = Line(center(l)\basis,  axis(l)\basis)
function rotate(l::Line, n::StaticVector{3}, θ::Real)
    # Rodrigues' rotation formula (rotate a line around a unit vector n by an angle θ)
    s, c = sincos(θ)
    cntr′::SVector{3,Float64} = c*center(l) + s*(n×center(l)) + (1-c)*(n⋅center(l)) * n
    axis′::SVector{3,Float64} = c*axis(l) + s*(n×axis(l)) + (1-c)*(n⋅axis(l)) * n
    return Line(cntr′, axis′)
end

# ---------------------------------------------------------------------------------------- #

struct Cylinder
    line::Line
    radius::Float64
end
Cylinder(cntr::AbstractVector{<:Real}, axis::AbstractVector{<:Real}, radius::Real) = Cylinder(Line(cntr, axis), radius)
line(c::Cylinder)   = c.line
center(c::Cylinder) = center(line(c))
axis(c::Cylinder)   = axis(line(c))
radius(c::Cylinder) = c.radius
cartesianize(c::Cylinder, basis::DirectBasis{3}) = Cylinder(transform(line(c), basis), radius(c))
latticize(c::Cylinder, basis::DirectBasis{3})    = Cylinder(transform(line(c), basis), radius(c))
rotate(c::Cylinder, n::StaticVector{3}, θ::Real) = Cylinder(rotate(line(c), n, θ), radius(c))

# ---------------------------------------------------------------------------------------- #

struct Sphere
    cntr::SVector{3,Float64}
    radius::Float64
end
center(s::Sphere) = s.cntr
radius(s::Sphere) = s.radius

# ---------------------------------------------------------------------------------------- #

function (*)(op::SymOperation{3}, c::Cylinder)
    cntr′ = rotation(op)*center(c) + translation(op)
    axis′ = rotation(op)*axis(c)
    axis′ = normalize(axis′)
    return Cylinder(Line(cntr′, axis′), radius(c))
end
function (==)(c1::Cylinder, c2::Cylinder)
    # axes are the same if they are colinear
    iszero(axis(c1) × axis(c2)) || return false
    # lines are the same if their centers lie on the line spanned by their axis
    iszero((center(c1) - center(c2)) × axis(c1)) || return false    
    # finally, cylinders are the same if they have the same radius
    return radius(c1) == radius(c2)
end
function isapprox(c1::Cylinder, c2::Cylinder;
            atol::Real=1e-12, rtol::Real=atol>0 ? 0.0 : √eps())
    norm(axis(c1) × axis(c2)) ≤ atol || return false
    norm((center(c1) - center(c2)) × axis(c1)) ≤ atol || return false  
    return isapprox(radius(c1), radius(c2); atol=atol, rtol=rtol)
end
Base.:+(c::Cylinder, v::StaticVector{3}) = Cylinder(Line(center(c)+v, axis(c)), radius(c))

# ---------------------------------------------------------------------------------------- #

function (*)(op::SymOperation{3}, s::Sphere)
    cntr′ = rotation(op)*center(s) + translation(op)
    return Sphere(cntr′, radius(s))
end
(==)(s1::Sphere, s2::Sphere) = center(s1) == center(s2) && radius(s1) == radius(s2)
function isapprox(s1::Sphere, s2::Sphere;
            atol::Real=1e-12, rtol::Real=atol>0 ? 0.0 : √eps())
    return isapprox(center(s1), center(s2); atol, rtol) && 
           isapprox(radius(s1), radius(s2); atol, rtol)
end
Base.:+(s::Sphere, v::StaticVector{3}) = Sphere(center(s)+v, radius(s))

# ---------------------------------------------------------------------------------------- #

function distance_to_cylinder_axis(r::AbstractVector{<:Real}, c::Cylinder)
    # cf. https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
    tmp = r - center(c)
    return norm(tmp - dot(tmp, axis(c))*axis(c))
end
function signed_distance(r::AbstractVector{<:Real}, c::Cylinder)
    # evaluates the shortest distance from an input point `r` to a cylinder `c`: if `r` is
    # outside `c`, the returned value is positive; conversely, if `r` is inside `c`, the
    # returned value is negative
    return distance_to_cylinder_axis(r, c) - radius(c)
end
function signed_distance(r::AbstractVector{<:Real}, cs::AbstractVector{<:Union{Cylinder, Sphere}})
    # evaluates the minimum signed distance from `r` to any number of cylinders or spheres `cs`
    d = Inf
    for c in cs
        d = min(d, signed_distance(r, c))
    end
    return d
end

function signed_distance(r::AbstractVector{<:Real}, s::Sphere)
    return norm(r - center(s)) - radius(s)
end

## --------------------------------------------------------------------------------------- #

include("cylinder_polygon_intersect.jl")

# ---------------------------------------------------------------------------------------- #

#=
function intersects(b::Box, s::Sphere)
    b_min, b_max = minimum(b), maximum(b)
    s_cntr = center(s)
    r = radius(s)
    for i in 1:3
        s_cntr[i] + r - b_min[i] > 0.0 || return false
        s_cntr[i] - r - b_max[i] < 0.0 || return false
    end
    return true
end
=#

# ---------------------------------------------------------------------------------------- #

function facets(uc :: Cell) # NB: we always return in _cartesian_ coordinates
    vs = Brillouin.vertices(uc)
    vs = if Brillouin.setting(uc) == Brillouin.LATTICE
        cartesianize.(vs, Ref(Brillouin.basis(uc))) # non-mutating, intentionally
    end
    return [vs[f] for f in Brillouin.faces(uc)]
end

function facets(Vs :: DirectBasis{3}; origin::SVector{3,Float64}=SVector(0.0,0.0,0.0))
    # returns a vector of vertices (i.e., faces, ordered counter-clockwise, i.e., with 
    # normal point "outward") of the unit cell's facets
    # - `origin`: where to place the center of the unit cell
    z = origin - sum(Vs) / 2
    A1 = [z, z + Vs[3], z + Vs[2] + Vs[3], z + Vs[2]]
    A2 = [z, z + Vs[1], z + Vs[3] + Vs[1], z + Vs[3]]
    A3 = [z, z + Vs[2], z + Vs[1] + Vs[2], z + Vs[1]]
    B1 = reverse(A1) .+ (Vs[1],)
    B2 = reverse(A2) .+ (Vs[2],)
    B3 = reverse(A3) .+ (Vs[3],)
    return [A1, A2, A3, B1, B2, B3]
end
# ---------------------------------------------------------------------------------------- #

const NEIGHBOR_TRANSLATIONS = 
    NTuple{26, SVector{3, Int}}([SVector(I[1],I[2],I[3]) for I in vec(CartesianIndices((-1:1, -1:1, -1:1))) if !iszero(I)])

function symmetrize(
    ops::AbstractVector{SymOperation{3}},
    cs::Union{Cylinder, Sphere, AbstractVector{Cylinder}, AbstractVector{Sphere}},
    boundary::Union{DirectBasis{3}, Cell{3}},
    fs::AbstractVector{<:AbstractVector{<:StaticVector{3, Float64}}} = facets(boundary);
    add_neighbors::Bool=true,
    cartesian_ops::Bool=false
)
    Rs = if boundary isa DirectBasis{3}
        boundary
    elseif boundary isa Cell{3}
        DirectBasis{3}(Brillouin.basis(boundary))
    else
        error("unreachable: `boundary` is not a `DirectBasis` or `Cell`")
    end

    cs isa AbstractVector || (cs = [cs])
    cs′ = copy(cs)
    for op in ops
        opᶜ = cartesian_ops ? op : cartesianize(op, Rs)
        for c in cs
            c′ = opᶜ * c
            isapproxin(c′, cs′) && continue # already seen
            any(f -> intersects(c′, f), fs) && push!(cs′, c′)
            
            add_neighbors || continue
            for t in NEIGHBOR_TRANSLATIONS
                c′′ = c′ + cartesianize(t, Rs)
                isapproxin(c′′, cs′) && continue # already seen
                any(f -> intersects(c′′, f), fs) && push!(cs′, c′′)
            end
        end
    end

    if length(cs′) == length(cs)
        # no new objects added ⇒ converged       
        # however, before returning we do the following check: in case the initial seed
        # point was outside the unit cell, there might be an element in `cs′` that really
        # doesn't intersect the unitcell; no point in keeping it around - filter it out here
        return filter!(c->any(f->intersects(c, f), fs), cs′) # no new cylinders added ⇒ converged
    else
        return symmetrize(ops, cs′, boundary, fs; add_neighbors=false, cartesian_ops)
    end
end

# ---------------------------------------------------------------------------------------- #
end # Module