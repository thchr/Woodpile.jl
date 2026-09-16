# geometric utilities shared by the primitive-polygon intersection tooling

using StaticArrays, LinearAlgebra

const INTERSECTION_DEFAULT_ATOL = 1e-10 # default absolute tolerance for intersection tests
norm2(v) = dot(v, v)

# calculate polygon normal using Newell's method, see e.g.:
# https://math.stackexchange.com/questions/2885839/ (the basic idea is to estimate the
# normal vector from the Stokes theorem, which relates the "area vector" A (which points
# along the normal vector!) according to A = ½ ∮ r·dr, with the integral looping around
# the polygon boundary (the area vector can be derived from Stokes theorem ∮F·dr=∫∇×F·dA,
# by picking F = c×r for some arbitrary constant vector c))
function face_normal(vs::AbstractVector{T}) where T <: StaticVector{3, <:Real}
    length(vs) < 3 && error("Polygon must have at least 3 vertices.")

    n = zero(SVector{3, Float64})
    N = length(vs)
    vᵢ₊₁ = vs[1]
    for i in eachindex(vs)
        vᵢ = vᵢ₊₁
        vᵢ₊₁ = vs[mod1(i+1, N)] # wraps around
        n += SVector{3, Float64}(
            (vᵢ[2] - vᵢ₊₁[2]) * (vᵢ[3] + vᵢ₊₁[3]),
            (vᵢ[3] - vᵢ₊₁[3]) * (vᵢ[1] + vᵢ₊₁[1]),
            (vᵢ[1] - vᵢ₊₁[1]) * (vᵢ[2] + vᵢ₊₁[2])
        )
    end

    A = norm(n) # this is (twice) the area of the polygon; `n` is 2× the area vector 𝐀=A𝐧
    iszero(A) && error("failed to find a consistent normal vector: polygon area is zero (vertices may be collinear)")
    return T(n / A)
end

# check if point `p` is inside a *convex* polygon (plane defined by vertices/normal)
function is_inside_convex_polygon(
    p::StaticVector{3,<:Real},
    vs::AbstractVector{<:StaticVector{3,<:Real}},
    n::StaticVector{3,<:Real} = face_normal(vs)
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

function _is_planar(
    vs::AbstractVector{<:StaticVector{3,<:Real}},
    n::StaticVector{3,<:Real} = face_normal(vs),
    atol::Real = INTERSECTION_DEFAULT_ATOL
    )
    v₁ = first(vs)
    for v in @views vs[2:end]
        δ = abs(dot(n, v - v₁))
        δ > atol && return false # point is not on the plane
    end
    return true
end
