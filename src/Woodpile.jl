module Woodpile

# ---------------------------------------------------------------------------------------- #

using Crystalline: SymOperation, rotation, translation, isapproxin, @S_str, compose
using StaticArrays
using LinearAlgebra

import Base: minimum, maximum, *, ==, isapprox

# ---------------------------------------------------------------------------------------- #

export Cylinder, Ray, Box, Sphere
export center, axis, radius
export symmetrize

# ---------------------------------------------------------------------------------------- #

struct Ray
    cntr::SVector{3,Float64}
    axis::SVector{3,Float64}
    function Ray(cntr::SVector{3,Float64}, axis::SVector{3,Float64})
        l = norm(axis)
        iszero(l) && error("Ray axis cannot be a zero vector")
        new(cntr, isone(l) ? axis : axis./l)
    end
end
Ray(cntr, axis) = Ray(convert(SVector{3,Float64}, cntr), convert(SVector{3,Float64}, axis))
center(r::Ray) = r.cntr
axis(r::Ray) = r.axis
# Rodrigues' rotation formula (rotate a ray around a unit vector n by an angle θ)
function rotate(r::Ray, n::StaticVector{3}, θ::Real)
    s, c = sincos(θ)
    cntr′::SVector{3,Float64} = c*r.cntr + s*(n×r.cntr) + (1-c)*(n⋅r.cntr)*n
    axis′::SVector{3,Float64} = c*r.axis + s*(n×r.axis) + (1-c)*(n⋅r.axis)*n
    return Ray(cntr′, axis′)
end

struct Cylinder
    ray::Ray
    radius::Float64
end
Cylinder(cntr::AbstractVector{<:Real}, axis::AbstractVector{<:Real}, radius::Real) = Cylinder(Ray(cntr, axis), radius)
ray(c::Cylinder)    = c.ray
center(c::Cylinder) = center(ray(c))
axis(c::Cylinder)   = axis(ray(c))
radius(c::Cylinder) = c.radius
rotate(c::Cylinder, n::StaticVector{3}, θ::Real) = Cylinder(rotate(ray(c), n, θ), radius(c))

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
Base.:+(c::Cylinder, v::StaticVector{3}) = Cylinder(Ray(center(c)+v, axis(c)), radius(c))

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

# ---------------------------------------------------------------------------------------- #

struct RayBoxIntersection
    bool::Bool
    tmin::Float64
    tmax::Float64
end
Base.Bool(i::RayBoxIntersection) = i.bool

function intersects(b::Box, r::Ray)
    # branchless ray-box intersection following https://tavianator.com/2011/ray_box.html
    b_min, b_max = minimum(b), maximum(b)
    r_cntr = center(r)
    r_axis = axis(r)
    inv_r_axis = 1 ./ r_axis

    tx1 = (b_min[1] - r_cntr[1]) * inv_r_axis[1]
    tx2 = (b_max[1] - r_cntr[1]) * inv_r_axis[1]
    tmin = min(tx1, tx2)
    tmax = max(tx1, tx2)

    ty1 = (b_min[2] - r_cntr[2]) * inv_r_axis[2]
    ty2 = (b_max[2] - r_cntr[2]) * inv_r_axis[2]
    tmin = max(tmin, min(ty1, ty2))
    tmax = min(tmax, max(ty1, ty2))
    
    tz1 = (b_min[3] - r_cntr[3]) * inv_r_axis[3]
    tz2 = (b_max[3] - r_cntr[3]) * inv_r_axis[3]   
    tmin = max(tmin, min(tz1, tz2))
    tmax = min(tmax, max(tz1, tz2))

    return RayBoxIntersection(tmax ≥ tmin, tmin, tmax)
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

struct SphereBoxIntersection
    bool::Bool
end
Base.Bool(i::SphereBoxIntersection) = i.bool

function intersects(b::Box, s::Sphere)
    b_min, b_max = minimum(b), maximum(b)
    s_cntr = center(s)
    r = radius(s)
    for i in 1:3
        s_cntr[i] + r - b_min[i] > 0.0 || return SphereBoxIntersection(false)
        s_cntr[i] - r - b_max[i] < 0.0 || return SphereBoxIntersection(false)
    end
    return SphereBoxIntersection(true)
end
might_intersects(b::Box, s::Sphere) = intersects(b, s)
might_intersects(s::Sphere, b::Box=Box()) = might_intersects(b, s)

# ---------------------------------------------------------------------------------------- #
const NEIGHBOR_TRANSLATIONS = 
    NTuple{26, SVector{3, Int}}([SVector(I[1],I[2],I[3]) for I in vec(CartesianIndices((-1:1, -1:1, -1:1))) if !iszero(I)])

function symmetrize(ops::AbstractVector{SymOperation{3}},
            cs::Union{AbstractVector{Cylinder}, AbstractVector{Sphere}};
            add_neighbors::Bool=true)

    cs′ = copy(cs)
    for op in ops
        for c in cs
            c′ = op * c
            (!isapproxin(c′, cs′) && Bool(might_intersects(c′))) && push!(cs′, c′)
            if add_neighbors
                for n in NEIGHBOR_TRANSLATIONS
                    c′′ = c′ + n
                    (!isapproxin(c′′, cs′) && Bool(might_intersects(c′′))) && push!(cs′, c′′)
                end
            end
        end
    end

    if length(cs′) == length(cs)
        # no extra objects were added, we converged       
        # however, before returning we do the following check: in case the initial seed
        # point was outside the unit cell, there might be an element in `cs′` that really
        # doesn't intersect the unitcell; no point in keeping it around - filter it out here
        return filter!(c->Bool(might_intersects(c)), cs′) # no extra cylinders were added, we converged
    else
        return symmetrize(ops, cs′; add_neighbors=false)
    end
end
symmetrize(ops::AbstractVector{SymOperation{3}}, c::Union{Cylinder, Sphere}; kwargs...) = symmetrize(ops, [c]; kwargs...)

function box_boundary(r, limits, boxtol=1e-8)
    x,y,z = r
    (x < limits[1]+boxtol || x > limits[2]-boxtol) && return 1e20
    (y < limits[3]+boxtol || y > limits[4]-boxtol) && return 1e20
    (z < limits[5]+boxtol || z > limits[6]-boxtol) && return 1e20
    return -1e20
end

# ---------------------------------------------------------------------------------------- #
end # Module