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
    NTuple{124, SVector{3, Int}}([SVector(I[1],I[2],I[3]) for I in vec(CartesianIndices((-2:2, -2:2, -2:2))) if !iszero(I)])

"""
    symmetrize(
        ops::AbstractVector{SymOperation{3}},
        cs::Union{Cylinder, AbstractVector{Cylinder}},
        boundary::Union{DirectBasis{3}, Cell{3}};
        add_neighbors::Bool=true,
        cartesian_ops::Bool=false) --> cs′ :: AbstractVector{<:Union{Cylinder, Sphere}}

Given a list of symmetry operations `ops` and a cylinder (or list of cylinders) `cs`, apply
every operation in `ops` to `cs` and aggregate the resulting distinct cylinders that
intersects the interior of `boundary`. This `boundary` can be specified either as a
`DirectBasis{3}` (from Bravais.jl/Crystalline.jl; then interpreted as a trapezoidal
boundary) or a `Cell{3}` (from Brillouin.jl; a Wigner-Seitz unit cell).

Note that the cylinders specified in `cs` must be specified in a Cartesian basis: i.e.,
their axis center, orientation, and radius all refer to a Cartesian coordinate system (The
same coordinate system as `boundary`). The symmetry operations `ops` can be specified in
a cartesian or a lattice basis (see `cartesian_ops` keyword argument).

## Keyword arguments
- `add_neighbors` (default: `true`): if `true`, the list of symmetry operations is augmented
  with translations to neighboring unit cells, up to 2nd nearest neighbor.
- `cartesian_ops` (default: `false`): if `true`, the symmetry operations are interpreted as
  specified in a Cartesian basis. If `false`, the symmetry operations are converted
  internally to a Cartesian basis using information from `boundary`.

## Examples
It is expected that Woodpile is used in conjunction with Crystalline.jl in order to access
the symmetry operations of space groups. The following example illustrates how to generate a
symmetric woodpile structures in space group 14 (P2₁/c):

```jl
julia> using Crystalline, Woodpile

julia> ops = spacegroup(14)
SpaceGroup{3} ⋕14 (P2₁/c) with 4 operations:
 1
 {2₀₁₀|0,½,½}
 -1
 {m₀₁₀|0,½,½}

julia> Rs = DirectBasis{3}([1.0, 0.0, 0.0], [0.0, 1.4, 0.0], [-0.45, 0.0, 0.5]) # monoclinic
DirectBasis{3} (monoclinic):
 [1.0, 0.0, 0.0]
 [0.0, 1.4, 0.0]
 [-0.45, 0.0, 0.5]

julia> cs = Cylinder([0, 0, 0], Rs[1], 0.15) # going through origo, along R₁, radius 0.15
Cylinder(Line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), 0.15)

julia> cs′ = symmetrize(ops, cs, Rs)
 Cylinder(Line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), 0.15)
 Cylinder(Line([-0.225, 0.7, 0.25], [-1.0, 0.0, 0.0]), 0.15)
 Cylinder(Line([-1.7750000000000001, -0.7, -0.25], [-1.0, 0.0, 0.0]), 0.15)
 Cylinder(Line([-1.7750000000000001, 0.7, -0.25], [-1.0, 0.0, 0.0]), 0.15)
 Cylinder(Line([-2.225, -0.7, 0.25], [-1.0, 0.0, 0.0]), 0.15)
```

Rather than the trapezoidal unit cell, we can also specify the boundary as the Wigner-Seitz
unit cell, using Brillouin.jl's `wignerseitz`:
```jl
julia> using Brillouin

julia> uc = wignerseitz(Rs);

julia> uc_cs′ = symmetrize(ops, cs, uc)
 Cylinder(Line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), 0.15)
 Cylinder(Line([-0.225, 0.7, 0.25], [-1.0, 0.0, 0.0]), 0.15)
 Cylinder(Line([-1.7750000000000001, -0.7, -0.25], [-1.0, 0.0, 0.0]), 0.15)
 Cylinder(Line([-1.7750000000000001, 0.7, -0.25], [-1.0, 0.0, 0.0]), 0.15)
 Cylinder(Line([-2.225, -0.7, 0.25], [-1.0, 0.0, 0.0]), 0.15)
 Cylinder(Line([1.55, 0.0, 0.5], [1.0, 0.0, 0.0]), 0.15)
 Cylinder(Line([-1.55, 0.0, -0.5], [-1.0, 0.0, 0.0]), 0.15)
```

In this case, more cylinders are needed: in general, the number of cylinders may differ for
trapezoidal and Wigner-Seitz unit cells (but the cylinders' enclosed volume is invariant).

## Visualization
The "symmetrized" cylinders can be visualized in the unit cell associated with `boundary`
using Makie.jl. For example, the following illustrates a symmetric woodpile structure in
space group 14 in its trapezoidal and Wigner Seitz unit cells, respectively:

```jl
julia> using GLMakie

julia> plot(cs′, Rs) # trapezoidal unit cell

julia> plot(uc_cs′, uc) # Wigner-Seitz unit cell
```

The number of samples used to resolve the cylinder isosurfaces can be controlled using
the keyword `samples` in the `plot` function. Similarly, the isosurfaces can be plotted
as `:merged` or `:individual` (default: `:merged`) via the `style` keyword argument.
"""
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