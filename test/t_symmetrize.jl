using Woodpile
using Woodpile: facets, intersects_cell, is_inside_cell
using Crystalline: spacegroup
using Bravais: DirectBasis
using Brillouin: wignerseitz, cartesianize, in_wignerseitz
using StaticArrays: SVector

# canonical, order-independent summary of a set of primitives (`Tuple`s, unlike `SVector`s,
# are `isless`-comparable and so can be `sort`ed)
sphere_keys(ss) = sort([(Tuple(center(s))..., radius(s)) for s in ss])

const Rs_cubic = DirectBasis{3}([1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0])
const Rs_mono  = DirectBasis{3}([1.0,0.0,0.0], [0.0,1.4,0.0], [-0.45,0.0,0.5])

@testset "unit cell facets & containment" begin
    # the facet normals of a trapezoidal cell must point out of the cell for
    # `is_inside_cell` to work; check this against a direct lattice-coordinate criterion
    for Rs in (Rs_cubic, Rs_mono)
        fs = facets(Rs)
        @test is_inside_cell(SVector(0.0, 0.0, 0.0), fs) # cell is centered at origo

        for _ in 1:250
            x = SVector{3,Float64}(3 .* (rand(3) .- 0.5)) # lattice coords ∈ [-1.5, 1.5]³
            rᶜ = SVector{3,Float64}(sum(x[i] * Rs[i] for i in 1:3))
            @test is_inside_cell(rᶜ, fs) == all(xᵢ -> -0.5 ≤ xᵢ ≤ 0.5, x)
        end
    end

    # the Wigner-Seitz facets from Brillouin.jl must be oriented consistently too
    for Rs in (Rs_cubic, Rs_mono)
        uc = wignerseitz(Rs)
        ucᶜ = cartesianize(uc)
        for fs in (facets(uc), facets(ucᶜ)) # must agree for either `setting(uc)`
            @test is_inside_cell(SVector(0.0, 0.0, 0.0), fs)
            for _ in 1:250
                x = SVector{3,Float64}(3 .* (rand(3) .- 0.5))
                rᶜ = SVector{3,Float64}(sum(x[i] * Rs[i] for i in 1:3))
                @test is_inside_cell(rᶜ, fs) == in_wignerseitz(rᶜ, ucᶜ)
            end
        end
    end
end # @testset "unit cell facets & containment"

@testset "intersects_cell" begin
    fs = facets(Rs_cubic) # cube spanning [-½, ½]³

    # a cylinder is infinite: it overlaps the cell iff it intersects one of its facets
    @test intersects_cell(Cylinder([0,0,0], [0,0,1], 0.1), fs) == true
    @test intersects_cell(Cylinder([2,0,0], [0,0,1], 0.1), fs) == false

    # a sphere is finite: it may be strictly interior and touch no facet at all
    @test intersects_cell(Sphere([0.0, 0.0, 0.0], 0.1), fs) == true
    @test intersects_cell(Sphere([0.45, 0.0, 0.0], 0.1), fs) == true # straddles a facet
    @test intersects_cell(Sphere([0.55, 0.0, 0.0], 0.1), fs) == true # just barely reaches
    @test intersects_cell(Sphere([0.7, 0.0, 0.0], 0.1), fs) == false # too far outside
    @test intersects_cell(Sphere([0.0, 0.0, 0.0], 10.0), fs) == true # cell inside sphere
end # @testset "intersects_cell"

@testset "symmetrize: cylinders" begin
    # regression check of the example from the `symmetrize` docstring
    ops = spacegroup(14)
    c = Cylinder([0, 0, 0], Rs_mono[1], 0.15)

    cs′ = symmetrize(ops, c, Rs_mono)
    @test cs′ isa Vector{Cylinder}
    @test length(cs′) == 5
    @test c ∈ cs′

    uc_cs′ = symmetrize(ops, c, wignerseitz(Rs_mono))
    @test length(uc_cs′) == 7
    @test c ∈ uc_cs′
end # @testset "symmetrize: cylinders"

@testset "symmetrize: spheres" begin
    ops = spacegroup(221) # Pm-3m

    # a lone sphere at the origin lies entirely inside the cell and intersects none of its
    # facets: it must nonetheless be retained. This case cannot arise for a `Cylinder`
    # (being infinite, it always exits through a facet), and is why `Sphere`s need the
    # `intersects_cell` check rather than a bare facet-intersection check
    s = Sphere([0.0, 0.0, 0.0], 0.2)
    ss′ = symmetrize(ops, s, Rs_cubic)
    @test ss′ isa Vector{Sphere}
    @test ss′ == [s]

    # Wyckoff position 8g of SG 221 sits at (x,x,x); for x = ¼ all 8 spheres of the orbit
    # lie strictly inside the cubic cell, and no periodic image reaches into it
    ss′ = symmetrize(ops, Sphere([0.25, 0.25, 0.25], 0.1), Rs_cubic)
    @test length(ss′) == 8
    @test all(s -> radius(s) == 0.1, ss′)
    @test sphere_keys(ss′) == sort([(x, y, z, 0.1) for x in (-0.25, 0.25)
                                                   for y in (-0.25, 0.25)
                                                   for z in (-0.25, 0.25)])

    # same, for a Wigner-Seitz boundary (identical to the cube, for a cubic lattice)
    ss′_ws = symmetrize(ops, Sphere([0.25, 0.25, 0.25], 0.1), wignerseitz(Rs_cubic))
    @test sphere_keys(ss′_ws) == sphere_keys(ss′)

    # a sphere at the cell corner (⅟₂,⅟₂,⅟₂) (Wyckoff 1b) has 8 periodic images reaching
    # into the cell, one at each corner
    ss′ = symmetrize(ops, Sphere([0.5, 0.5, 0.5], 0.1), Rs_cubic)
    @test length(ss′) == 8
    @test all(s -> all(x -> abs(x) ≈ 0.5, center(s)), ss′)

    # spheres in a lower-symmetry, non-orthogonal setting: the orbit of a general position
    # of SG 14 has 4 members, all lying inside the cell (for a small enough radius)
    fs_mono = facets(Rs_mono)
    ss′ = symmetrize(spacegroup(14), Sphere([0.1, 0.2, 0.05], 0.04), Rs_mono)
    @test length(ss′) == 4
    @test allunique(ss′)
    @test all(s -> is_inside_cell(center(s), fs_mono), ss′)

    # at a larger radius, two periodic images of that orbit poke into the cell through the
    # R₃ facets: they are centered outside the cell, but must be retained nonetheless
    ss′ = symmetrize(spacegroup(14), Sphere([0.1, 0.2, 0.05], 0.06), Rs_mono)
    @test length(ss′) == 6
    @test count(s -> !is_inside_cell(center(s), fs_mono), ss′) == 2
end # @testset "symmetrize: spheres"

@testset "symmetrize: mixed cylinders & spheres" begin
    ops = spacegroup(221)
    c = Cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.1)
    s = Sphere([0.25, 0.25, 0.25], 0.1)

    cs′ = symmetrize(ops, Primitive[c, s], Rs_cubic)
    @test cs′ isa Vector{Primitive}
    @test sphere_keys(filter(x -> x isa Sphere, cs′)) ==
          sphere_keys(symmetrize(ops, s, Rs_cubic))
    @test length(filter(x -> x isa Cylinder, cs′)) ==
          length(symmetrize(ops, c, Rs_cubic))
end # @testset "symmetrize: mixed cylinders & spheres"
