using Woodpile
using Woodpile: Line, line
using Bravais: DirectBasis, cartesianize, latticize
using LinearAlgebra: normalize
using StaticArrays: SVector

@testset "cartesianize & latticize" begin
    Rs = DirectBasis{3}([1.0,0.0,0.0], [0.0,1.4,0.0], [-0.45,0.0,0.5]) # monoclinic
    Is = DirectBasis{3}([1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0])   # Cartesian

    # a `Line`, with center & axis referred to the lattice basis
    l = Line([0.25, 0.5, 0.25], [0.0, 0.0, 1.0])
    lᶜ = cartesianize(l, Rs)
    @test center(lᶜ) ≈ cartesianize(center(l), Rs)
    @test axis(lᶜ) ≈ normalize(cartesianize(axis(l), Rs)) # a `Line` normalizes its axis
    @test center(latticize(lᶜ, Rs)) ≈ center(l) # round-trip
    @test axis(latticize(lᶜ, Rs)) ≈ axis(l)

    c = Cylinder([0.25, 0.5, 0.25], [0.0, 0.0, 1.0], 0.15)
    cᶜ = cartesianize(c, Rs)
    @test center(cᶜ) ≈ center(cartesianize(line(c), Rs))
    @test radius(cᶜ) == radius(c) # radii are lengths, and are left untouched
    @test latticize(cᶜ, Rs) ≈ c   # round-trip

    s = Sphere([0.25, 0.5, 0.25], 0.15)
    sᶜ = cartesianize(s, Rs)
    @test center(sᶜ) ≈ cartesianize(center(s), Rs)
    @test center(sᶜ) ≈ SVector(0.25 - 0.45*0.25, 1.4*0.5, 0.5*0.25)
    @test radius(sᶜ) == radius(s)
    @test latticize(sᶜ, Rs) ≈ s # round-trip

    # a Cartesian basis is a fixed point of both
    @test cartesianize(c, Is) ≈ c
    @test latticize(c, Is) ≈ c
    @test cartesianize(s, Is) ≈ s
    @test latticize(s, Is) ≈ s
end # @testset "cartesianize & latticize"
