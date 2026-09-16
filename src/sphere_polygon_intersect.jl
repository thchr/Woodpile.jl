# minimum squared distance between a point `p` and the segment from `A` to `B`: the closest
# point on the segment is Q(t) = A + t(B-A) with t = dot(p-A, B-A)/|B-A|², clamped to [0, 1]
function segment2point_dist2(
    A::StaticVector{3,<:Real},
    B::StaticVector{3,<:Real},
    p::StaticVector{3,<:Real},
    atol::Real = INTERSECTION_DEFAULT_ATOL
)
    Δ = B - A
    c = norm2(Δ)
    c < atol && return norm2(p - A) # degenerate segment (a point)
    t′ = clamp(dot(p - A, Δ) / c, 0.0, 1.0)
    return norm2(p - (A + t′*Δ))
end

# minimum squared distance between a point `p` and the boundary (i.e., the edges) of a
# polygon, with early exit if any distance is less than R² + atol (i.e., intersection occurs);
# NB: the polygon's vertices need not be checked separately, since they are edge end-points
function polygon2point_dist2(
    vs::AbstractVector{<:StaticVector{3,<:Real}},
    p::StaticVector{3,<:Real},
    R²::Real,
    atol::Real = INTERSECTION_DEFAULT_ATOL
)
    min_d² = typemax(Float64) # initialize to max value
    N = length(vs)
    vᵢ₊₁ = vs[1]
    for i in eachindex(vs)
        vᵢ = vᵢ₊₁
        vᵢ₊₁ = vs[mod1(i + 1, N)] # wrap-around
        d² = segment2point_dist2(vᵢ, vᵢ₊₁, p, atol)
        d² ≤ R² + atol && return d² # early exit (w/ tolerance for comparison)
        min_d² = min(min_d², d²)
    end

    return min_d²
end

## --------------------------------------------------------------------------------------- #

"""
    intersects(s::Sphere, vs::AbstractVector{<:StaticVector{3,<:Real}};
               atol::Real = INTERSECTION_DEFAULT_ATOL)

Determine if a `Sphere`, `s`, intersects a flat *convex* polygon, as specified by its
(clockwise or counter-clockwise ordered) vertices `vs`.

Returns `true` if an intersection exists and `false` otherwise. Errors on invalid input.

## Optional arguments
- `atol`: tolerance for floating-point comparisons.
"""
function intersects(
    s::Sphere,
    vs::AbstractVector{<:StaticVector{3,<:Real}};
    atol::Real = INTERSECTION_DEFAULT_ATOL
)
    length(vs) < 3 && error("polygon must have at least 3 vertices.")

    # sphere properties
    p = center(s)
    R = radius(s)
    R² = R * R

    # polygon plane properties
    n = face_normal(vs)
    v₁ = first(vs)
    _is_planar(vs, n, atol) || error("polygon vertices are not coplanar within tolerance atol.")

    d_plane = dot(n, p - v₁) # signed distance from sphere center to the polygon's plane
    abs(d_plane) > R + atol && return false # sphere doesn't reach the plane of the polygon

    # the sphere reaches the plane: its closest approach to the polygon is either at the
    # projection of its center onto the plane (if that projection is inside the polygon) or
    # at the polygon's boundary (if it is not)
    p′ = p - d_plane * n # projection of sphere center onto the polygon's plane
    is_inside_convex_polygon(p′, vs, n) && return true

    return polygon2point_dist2(vs, p, R², atol) ≤ R² + atol
end

## --------------------------------------------------------------------------------------- #

function intersects_cell(
    s::Sphere,
    fs::AbstractVector{<:AbstractVector{<:StaticVector{3,<:Real}}};
    atol::Real = INTERSECTION_DEFAULT_ATOL
)
    return is_inside_cell(center(s), fs; atol) || any(f -> intersects(s, f; atol), fs)
end
