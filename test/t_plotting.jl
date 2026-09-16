using Woodpile
using Crystalline: spacegroup
using Bravais: DirectBasis
using Brillouin: wignerseitz
using Makie: Makie, FigureAxisPlot, GeometryBasics

# per-dimension extent of the isosurface meshed by `plot(...)`: used below as a cheap proxy
# for "the primitives were actually meshed, and at roughly the right place"
function mesh_extrema(p)
    vs = GeometryBasics.coordinates(Makie.to_value(p[1]))
    isempty(vs) && return nothing
    return (extrema(v->v[1], vs), extrema(v->v[2], vs), extrema(v->v[3], vs))
end

@testset "plotting (Makie extension)" begin
    Rs = DirectBasis{3}([1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0])
    ops = spacegroup(221)

    c = Cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.15)
    s = Sphere([0.25, 0.25, 0.25], 0.15)
    cs′ = symmetrize(ops, c, Rs)
    ss′ = symmetrize(ops, s, Rs)
    @test !isempty(ss′)

    @testset "$(nameof(typeof(boundary))) boundary" for boundary in (Rs, wignerseitz(Rs))
        for style in (:merged, :individual)
            # spheres alone
            fap = Makie.plot(ss′, boundary; samples=40, style=style)
            @test fap isa FigureAxisPlot
            ex = mesh_extrema(fap.plot)
            @test ex !== nothing # a non-empty isosurface was meshed
            if style == :merged # (`:individual` returns only the *last* mesh plot)
                # the 8 spheres of radius 0.15 centered at (±¼,±¼,±¼) span [-0.4, 0.4] in
                # each direction; allow a margin for the finite marching-cubes resolution
                for (lo, hi) in ex
                    @test -0.45 ≤ lo ≤ -0.3
                    @test  0.3  ≤ hi ≤  0.45
                end
            end

            # cylinders alone (regression: plotting must be unaffected by `Sphere` support)
            fap = Makie.plot(cs′, boundary; samples=40, style=style)
            @test mesh_extrema(fap.plot) !== nothing

            # cylinders and spheres together
            fap = Makie.plot(Primitive[cs′..., ss′...], boundary; samples=40, style=style)
            @test mesh_extrema(fap.plot) !== nothing
        end

        @test_throws "unsupported style" Makie.plot(ss′, boundary; samples=10, style=:nope)
    end

    # a boundary-free plot (uses a unit cube in lattice coordinates)
    fap = Makie.plot(ss′; samples=40)
    @test fap isa FigureAxisPlot
    @test mesh_extrema(fap.plot) !== nothing
end # @testset "plotting (Makie extension)"
