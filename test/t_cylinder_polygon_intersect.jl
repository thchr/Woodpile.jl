using Woodpile: Cylinder, intersects

@testset "cylinder-polygon intersection" begin

    c12 = Cylinder([0,0,5], [1,0,0], 0.5) # axis along x, a z-height 5, radius 0.5
    poly1 = SVector{3, Float64}[ # square on the xy-plane centered at (2,0,0)
        [1.5, -0.5, 0.0], [2.5, -0.5, 0.0], [2.5, 0.5, 0.0], [1.5, 0.5, 0.0]]
    @test intersects(c12, poly1) == false

    poly2 = poly1 .- Ref(SVector(2, 0.0, -4.75)) # like `poly1` but centered at (0,0,4.75)
    @test intersects(c12, poly2) == true

    c3 = Cylinder([1.2, 0, 0], [0,0,1], 0.5) # axis along z, at x-height 1.2, radius 0.5
    square = SVector{3, Float64}[ # square at xy-plane, centered at origo
        [-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [-1.0, 1.0, 0.0]]
    @test intersects(c3, square) == true # (proximity check)

    # axis parallel with polygon plane and close enough to intersect
    c4 = Cylinder([0, 0, 0.4], [1, 0, 0], 0.5) # axis along x, at z-height 0.4, radius 0.5
    @test intersects(c4, square) == true # (parallel & proximity check passes)

    # axis parallel with polygon plane and too far away to intersect
    c5 = Cylinder([0, 0, 1.0], [1, 0, 0], 0.5) # axis along x, at z-height 1, radius 0.5
    @test intersects(c5, square) == false # (parallel & d_plane > R)

    # axis intersects polygon
    @test intersects(Cylinder([0,0,0], [0,0,1], 0.5), square) == true

    # planarity violated
    poly6 = SVector{3, Float64}[[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0], 
                                [-1.0, 1.0, 0.1] #= z != 0 =#] 
    @test_throws "vertices are not coplanar" intersects(c3, poly6)

    # collinear vertices in polygon
    poly7 = SVector{3, Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
    @test_throws "failed to find a consistent normal vector" intersects(c3, poly7)

end # @testset "cylinder-polygon intersection"