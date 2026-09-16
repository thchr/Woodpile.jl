using Woodpile: Sphere, intersects
using StaticArrays: SVector

square = SVector{3, Float64}[ # side-2 square in the xy-plane, centered at origo
    [-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [-1.0, 1.0, 0.0]]

@testset "sphere-polygon intersection" begin

    # sphere pierced by the polygon's interior
    @test intersects(Sphere([0.0, 0.0, 0.0], 0.5), square) == true
    @test intersects(Sphere([0.5, -0.5, 0.25], 0.5), square) == true

    # sphere too far from the polygon's plane (projection would be inside the polygon)
    @test intersects(Sphere([0.0, 0.0, 1.5], 0.5), square) == false
    @test intersects(Sphere([0.0, 0.0, -1.5], 0.5), square) == false

    # sphere close enough to the plane, but its projection falls outside the polygon and it
    # is too far from the polygon's boundary
    @test intersects(Sphere([2.0, 0.0, 0.25], 0.5), square) == false
    @test intersects(Sphere([2.0, 2.0, 0.0], 0.5), square) == false # (outside a corner)

    # sphere whose projection falls outside the polygon, but which reaches an edge
    @test intersects(Sphere([1.3, 0.0, 0.25], 0.5), square) == true
    # ... or a corner (at distance √2·0.3 ≈ 0.424 from the corner (1,1,0))
    @test intersects(Sphere([1.3, 1.3, 0.0], 0.5), square) == true
    @test intersects(Sphere([1.3, 1.3, 0.0], 0.4), square) == false

    # tangency: sphere touching the polygon's plane, interior, edge, and corner
    @test intersects(Sphere([0.0, 0.0, 0.5], 0.5), square) == true
    @test intersects(Sphere([1.5, 0.0, 0.0], 0.5), square) == true  # touches edge x=1
    @test intersects(Sphere([1.5, 0.0, 0.0], 0.49), square) == false
    @test intersects(Sphere([1.0, 1.0, 0.5], 0.5), square) == true  # touches corner (1,1,0)
    @test intersects(Sphere([1.0, 1.0, 0.5], 0.49), square) == false

    # `intersects` is a *solid* ball-vs-polygon test: a sphere that swallows the polygon
    # entirely intersects it, even though the sphere's *surface* does not
    @test intersects(Sphere([0.0, 0.0, 0.0], 100.0), square) == true
    @test intersects(Sphere([0.0, 0.0, 50.0], 100.0), square) == true

    # triangles work too (i.e., we don't assume 4 vertices)
    triangle = SVector{3, Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    @test intersects(Sphere([0.2, 0.2, 0.1], 0.2), triangle) == true
    @test intersects(Sphere([1.0, 1.0, 0.0], 0.2), triangle) == false
    @test intersects(Sphere([1.0, 1.0, 0.0], 0.8), triangle) == true # reaches hypotenuse

    # error paths (shared with the cylinder method)
    poly6 = SVector{3, Float64}[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0], 
                                [-1.0, 1.0, 0.1] #= z != 0 =#]
    @test_throws "vertices are not coplanar" intersects(Sphere([0.0,0.0,0.0], 0.5), poly6)
    poly7 = SVector{3, Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
    @test_throws "failed to find a consistent normal vector" intersects(Sphere([0.0,0.0,0.0], 0.5), poly7)
    @test_throws "at least 3 vertices" intersects(Sphere([0.0,0.0,0.0], 0.5), square[1:2])

end # @testset "sphere-polygon intersection"
