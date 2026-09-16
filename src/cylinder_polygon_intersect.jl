# squared distance between point p and line defined by `l`
function point2line_dist2(p::StaticVector{3,Float64}, l::Line)
    w = p - center(l)
    d² = norm2(w) - dot(w, axis(l))^2
    return max(0.0, d²) # clamp for numerical robustness
end

# Minimum squared distance between segment from `A` to `B` and line `l`:
# The squared distance between the axis P(s) = center(l) + s·axis(l) and the line associated 
# with the segment, Q(t) = `A` + t·(B-A), is a quadratic function in t and s. 
# The minimum of this quadratic function is found at some values (t′, s′), which can be 
# computed from Δ = B-A and w = A-center(l), and b = dot(axis(l), v), c = dot(Δ, Δ),
# d = dot(axis, w), and e = dot(Δ, w), constants that all appear in the quadratic form.
# Finally, to get the _segment_ minimum distance, we clamp t′ to the range [0, 1].
function segment2line_dist2(
    A::StaticVector{3,Float64},
    B::StaticVector{3,Float64},
    l::Line,
    atol::Float64 = INTERSECTION_DEFAULT_ATOL
)
    Δ = B - A
    w = A - center(l)

    # a = dot(axis(l), axis(l)) = 1.0 (normalized)
    b = dot(axis(l), Δ)
    c = norm2(Δ)

    # check if segment is degenerate (a point)
    c < atol && return point2line_dist2(A, l)

    d = dot(axis(l), w)
    e = dot(Δ, w)
    _det = c - b*b # = a*c - b*b = c - b*b since a=1
    if abs(_det) < atol # segment direction is parallel to line direction
        # distance is constant along the line containing the segment
        return point2line_dist2(A, l)
    end

    # calculate segment parameter 't' for closest approach point: derived by minimizing
    # ``|(A + t*Δ) - (center(l) + s*axis(l))|²``, leading to `t = (d*b - e)/(c - b²)`
    t = (b*e - c*d) / _det # for entire line (not necessarily segment)
    t = (d*b - e) / _det
    t′ = clamp(t, 0.0, 1.0) # clamp `t` to [0, 1] for the segment
    closest_point_on_segment = A + t′ * Δ

    # distance from this segment point to the infinite line (`axis(l)`)
    return point2line_dist2(closest_point_on_segment, l)
end

# minimum squared distance between polygon vertices/edges and line (`axis(l)`), with early
# exit if any distance is less than R² + atol (i.e., intersection occurs)
function polygon2line_dist2(
    vs::AbstractVector{<:StaticVector{3,Float64}},
    l::Line,
    R²::Float64,
    atol::Float64 = INTERSECTION_DEFAULT_ATOL
)    
    # check vertices
    min_d² = typemax(Float64) # initialize to max value
    for v in vs
        d² = point2line_dist2(v, l)
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
        d² = segment2line_dist2(vᵢ, vᵢ₊₁, l, atol)
        d² ≤ R² + atol && return d² # early exit (w/ tolerance for comparison)
        min_d² = min(min_d², d²)
    end

    return min_d²
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
    l = line(c)
    R = radius(c)
    R² = R * R

    # polygon plane properties
    n = face_normal(vs)
    v₁ = first(vs)
    _is_planar(vs, n, atol) || error("polygon vertices are not coplanar within tolerance atol.")

    dot_axis_normal = dot(axis(l), n)
    if abs(dot_axis_normal) < atol
        # branch 1: axis is parallel to the plane of the face
        d_plane = abs(dot(n, center(l) - v₁)) # distance from cylinder axis to plane

        if d_plane > R + atol
            return false # cylinder axis too far from plane
        else
            # check proximity of polygon vertices/edges to cylinder axis
            return polygon2line_dist2(vs, l, R², atol) ≤ R² + atol
        end

    else
        # branch 2: axis intersects the plane of the face; find intersection point
        t = dot(n, v₁ - center(l)) / dot_axis_normal
        P_intersect = center(l) + t * axis(l)

        # check if the intersection point is inside the convex polygon
        if is_inside_convex_polygon(P_intersect, vs, n)
            return true # axis pierces the polygon face or boundary
        else
            # axis intersects plane outside polygon; check proximity
            return polygon2line_dist2(vs, l, R², atol) ≤ R² + atol
        end
    end
end

## --------------------------------------------------------------------------------------- #

function intersects_cell(
    c::Cylinder,
    fs::AbstractVector{<:AbstractVector{<:StaticVector{3,Float64}}};
    atol::Float64 = INTERSECTION_DEFAULT_ATOL
)
    return any(f -> intersects(c, f; atol), fs)
end
