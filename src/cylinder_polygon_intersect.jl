using StaticArrays, LinearAlgebra

const INTERSECTION_DEFAULT_ATOL = 1e-10 # default absolute tolerance for intersection tests
norm2(v) = dot(v, v)

# calculate polygon normal using Newell's method, see e.g.:
# https://math.stackexchange.com/questions/2885839/ (the basic idea is to estimate the
# normal vector from the Stokes theorem, which relates the "area vector" A (which points
# along the normal vector!) according to A = ½ ∮ r·dr, with the integral looping around
# the polygon boundary (the area vector can be derived from Stokes theorem ∮F·dr=∫∇×F·dA,
# by picking F = c×r for some arbitrary constant vector c))
function face_normal(vs::AbstractVector{T}) where T <: StaticVector{3,Float64}
    length(vs) < 3 && error("Polygon must have at least 3 vertices.")

    n = zero(SVector{3, Float64})
    N = length(vs)
    vᵢ₊₁ = vs[1]
    for i in eachindex(vs)
        vᵢ = vᵢ₊₁
        vᵢ₊₁ = vs[mod1(i+1, N)] # wraps around
        n += SVector{3,Float64}(
            (vᵢ[2] - vᵢ₊₁[2]) * (vᵢ[3] + vᵢ₊₁[3]),
            (vᵢ[3] - vᵢ₊₁[3]) * (vᵢ[1] + vᵢ₊₁[1]),
            (vᵢ[1] - vᵢ₊₁[1]) * (vᵢ[2] + vᵢ₊₁[2])
        )
    end

    A = norm(n) # this is (twice) the area of the polygon; `n` is 2× the area vector 𝐀=A𝐧
    iszero(A) && error("failed to find a consistent normal vector: polygon area is zero (vertices may be collinear)")
    return T(n / A)
end

# squared distance between point p and line defined by ray
function point2ray_dist2(p::StaticVector{3,Float64}, r::Ray)
    w = p - center(r)
    d² = norm2(w) - dot(w, axis(r))^2
    return max(0.0, d²) # clamp for numerical robustness
end

# Minimum squared distance between segment from `A` to `B` and ray `r`:
# The squared distance between the axis P(s) = center(r) + s·axis(r) and the line associated 
# with the segment, Q(t) = `A` + t·(B-A), is a quadratic function in t and s. 
# The minimum of this quadratic function is found at some values (t′, s′), which can be 
# computed from Δ = B-A and w = A-center(r), and b = dot(axis(r), v), c = dot(Δ, Δ),
# d = dot(axis, w), and e = dot(Δ, w), constants that all appear in the quadratic form.
# Finally, to get the _segment_ minimum distance, we clamp t′ to the range [0, 1].
function segment2ray_dist2(
    A::StaticVector{3,Float64},
    B::StaticVector{3,Float64},
    r::Ray,
    atol::Float64 = INTERSECTION_DEFAULT_ATOL
)
    Δ = B - A
    w = A - center(r)

    # a = dot(axis(r), axis(r)) = 1.0 (normalized)
    b = dot(axis(r), Δ)
    c = norm2(Δ)

    # check if segment is degenerate (a point)
    c < atol && return point2ray_dist2(A, r)

    d = dot(axis(r), w)
    e = dot(Δ, w)
    _det = c - b*b # = a*c - b*b = c - b*b since a=1
    if abs(_det) < atol # segment direction is parallel to line direction
        # distance is constant along the line containing the segment
        return point2ray_dist2(A, r)
    end

    # calculate segment parameter 't' for closest approach point: derived by minimizing
    # ``|(A + t*Δ) - (center(r) + s*axis(r))|²``, leading to `t = (d*b - e)/(c - b²)`
    t = (b*e - c*d) / _det # for entire line (not necessarily segment)
    t = (d*b - e) / _det
    t′ = clamp(t, 0.0, 1.0) # clamp `t` to [0, 1] for the segment
    closest_point_on_segment = A + t′ * Δ

    # distance from this segment point to the infinite line (`axis(r)`)
    return point2ray_dist2(closest_point_on_segment, r)
end

# minimum squared distance between polygon vertices/edges and line (`axis(r)`), with early
# exit if any distance is less than R² + atol (i.e., intersection occurs)
function polygon2ray_dist2(
    vs::AbstractVector{<:StaticVector{3,Float64}},
    r::Ray,
    R²::Float64,
    atol::Float64 = INTERSECTION_DEFAULT_ATOL
)    
    # check vertices
    min_d² = typemax(Float64) # initialize to max value
    for v in vs
        d² = point2ray_dist2(v, r)
        if d² <= R² + atol # add tolerance for comparison
            return d² # early exit
        end
        d² ≤ R² + atol && return d² # early exit (w/ tolerance for comparison)
        min_d² = min(min_d², d²)
    end

    # check edges
    N = length(vs)
    vᵢ₊₁ = vs[1]
    for i in eachindex(vs)
        vᵢ = vᵢ₊₁
        vᵢ₊₁ = vs[mod1(i + 1, N)] # wrap-around
        d² = segment2ray_dist2(vᵢ, vᵢ₊₁, r, atol)
        d² ≤ R² + atol && return d² # early exit (w/ tolerance for comparison)
        min_d² = min(min_d², d²)
    end

    return min_d²
end

# check if point `p` is inside a *convex* polygon (plane defined by vertices/normal)
function is_inside_convex_polygon(
    p::StaticVector{3,Float64},
    vs::AbstractVector{<:StaticVector{3,Float64}},
    n::StaticVector{3,Float64} = face_normal(vs)
)
    # assumes vertices are ordered consistently (counter-clockwise ordering relative to `n`)
    # assumes point lies on the polygon's plane (checked before calling)
    N = length(vs)
    vᵢ₊₁ = vs[1]
    for i in eachindex(vs)
        vᵢ = vᵢ₊₁
        vᵢ₊₁ = vs[mod1(i + 1, N)] # wrap-around

        edge = vᵢ₊₁ - vᵢ
        Δᵢ = p - vᵢ # vector from edge start to point

        # edge normal pointing inwards (assuming counter-clockwise vertices relative to `n`)
        edge_n = cross(edge, n) # unnormalized edge normal; don't need correct scale

        # point must be on the negative side of the plane defined by edge & normal
        # (if dot-product is positive, `p` is outside the half-space defined by the `edge`)
        dot(Δᵢ, edge_n) > 0.0 && return false # point is outside this edge's half-space
    end

    # if inside or on boundary of all edge half-spaces, it's inside the convex polygon
    return true
end

## --------------------------------------------------------------------------------------- #

"""
    intersects(c::Cylinder, vs::AbstractVector{<:StaticVector{3,Float64}};
               atol::Float64 = INTERSECTION_DEFAULT_ATOL)

Determine if an infinitely extended `Cylinder`, `c`, intersects a flat *convex* polygon, as
specified by its (clockwise or counter-clockwise ordered) vertices `vs`.

Returns `true` if an intersection exists and `false` otherwise. Errors on invalid input.

## Optional arguments
- `atol`: tolerance for floating-point comparisons.
"""
function intersects(
    c::Cylinder,
    vs::AbstractVector{<:StaticVector{3,Float64}};
    atol::Float64 = INTERSECTION_DEFAULT_ATOL
)
    length(vs) < 3 && error("polygon must have at least 3 vertices.")

    # cylinder properties
    r = ray(c)
    R = radius(c)
    R² = R * R

    # polygon plane properties
    n = face_normal(vs)
    v₁ = first(vs)
    _is_planar(vs, n, atol) || error("polygon vertices are not coplanar within tolerance atol.")

    dot_axis_normal = dot(axis(r), n)
    if abs(dot_axis_normal) < atol
        # branch 1: axis is parallel to the plane of the face
        d_plane = abs(dot(n, center(r) - v₁)) # distance from cylinder axis to plane

        if d_plane > R + atol
            return false # cylinder axis too far from plane
        else
            # check proximity of polygon vertices/edges to cylinder axis
            return polygon2ray_dist2(vs, r, R², atol) ≤ R² + atol
        end

    else
        # branch 2: axis intersects the plane of the face; find intersection point
        t = dot(n, v₁ - center(r)) / dot_axis_normal
        P_intersect = center(r) + t * axis(r)

        # check if the intersection point is inside the convex polygon
        if is_inside_convex_polygon(P_intersect, vs, n)
            return true # axis pierces the polygon face or boundary
        else
            # axis intersects plane outside polygon; check proximity
            return polygon2ray_dist2(vs, r, R², atol) ≤ R² + atol
        end
    end
end

function _is_planar(
    vs::AbstractVector{<:StaticVector{3,Float64}},
    n::StaticVector{3,Float64} = face_normal(vs),
    atol::Float64 = INTERSECTION_DEFAULT_ATOL
    )
    v₁ = first(vs)
    for v in @views vs[2:end]
        δ = abs(dot(n, v - v₁))
        δ > atol && return false # point is not on the plane
    end
    return true
end