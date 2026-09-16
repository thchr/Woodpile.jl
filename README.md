# Woodpile.jl

Julia package to generate arbitrary woodpile-like structures with space group symmetry.

The functionality is provided via `symmetrize` (see `help> symmetrize`), which "symmetrizes" (i.e., adds symmetrically related) `Cylinder` and `Sphere` objects over a set of symmetry operations, retaining those that are inequivalent and inside the unit cell.

## Installation
Woodpile.jl is not registered in the Julia registry. To install, add via the GitHub URL:

```
pkg> add https://github.com/thchr/Woodpile.jl
```

Woodpile depends on and is expected to be used in combination with [Crystalline.jl](https://github.com/thchr/Crystalline.jl) and [Brillouin.jl](https://github.com/thchr/Brillouin.jl).
## Functionality
The docstring of `symmetrize` is copied below:

> ```jl
>     symmetrize(
>         ops::AbstractVector{SymOperation{3}},
>         cs::Union{Primitive, AbstractVector{<:Primitive}},
>         boundary::Union{DirectBasis{3}, Cell{3}};
>         add_neighbors::Bool=true,
>         cartesian_ops::Bool=false)
>     --> cs′ :: Vector{<:Primitive}
> ```
> Given a list of symmetry operations `ops` and a geometric primitive (or list of primitives)
> `cs`, apply every operation in `ops` to `cs` and aggregate the resulting distinct primitives
> that overlap the interior of `boundary`. This `boundary` can be specified either as a
> `DirectBasis{3}` (from Bravais.jl/Crystalline.jl; then interpreted as a trapezoidal
> boundary) or a `Cell{3}` (from Brillouin.jl; a Wigner-Seitz unit cell).
> 
> The primitives `cs` can be `Cylinder`s (infinitely extended) or `Sphere`s, and the two may
> be mixed in a single call (by giving `cs` as a `Vector{Primitive}`).
> 
> Note that the primitives specified in `cs` must be specified in a Cartesian basis: i.e.,
> their center, orientation, and radius all refer to a Cartesian coordinate system (the
> same coordinate system as `boundary`). The symmetry operations `ops` can be specified in
> a cartesian or a lattice basis (see `cartesian_ops` keyword argument).
> 
> ## Keyword arguments
> - `add_neighbors` (default: `true`): if `true`, the list of symmetry operations is augmented
>   with translations to neighboring unit cells, up to 2nd nearest neighbor.
> - `cartesian_ops` (default: `false`): if `true`, the symmetry operations are interpreted as
>   specified in a Cartesian basis. If `false`, the symmetry operations are converted
>   internally to a Cartesian basis using information from `boundary`.
> 
> ## Examples
> It is expected that Woodpile is used in conjunction with Crystalline.jl in order to access
> the symmetry operations of space groups. The following example illustrates how to generate a
> symmetric woodpile structures in space group 14 (P2₁/c):
> 
> ```jl
> julia> using Crystalline, Woodpile
> 
> julia> ops = spacegroup(14)
> SpaceGroup{3} ⋕14 (P2₁/c) with 4 operations:
>  1
>  {2₀₁₀|0,½,½}
>  -1
>  {m₀₁₀|0,½,½}
> 
> julia> Rs = DirectBasis{3}([1.0, 0.0, 0.0], [0.0, 1.4, 0.0], [-0.45, 0.0, 0.5]) # monoclinic
> DirectBasis{3} (monoclinic):
>  [1.0, 0.0, 0.0]
>  [0.0, 1.4, 0.0]
>  [-0.45, 0.0, 0.5]
> 
> julia> cs = Cylinder([0, 0, 0], Rs[1], 0.15) # going through origo, along R₁, radius 0.15
> Cylinder(Line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), 0.15)
> 
> julia> cs′ = symmetrize(ops, cs, Rs)
>  Cylinder(Line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), 0.15)
>  Cylinder(Line([-0.225, 0.7, 0.25], [-1.0, 0.0, 0.0]), 0.15)
>  Cylinder(Line([-1.7750000000000001, -0.7, -0.25], [-1.0, 0.0, 0.0]), 0.15)
>  Cylinder(Line([-1.7750000000000001, 0.7, -0.25], [-1.0, 0.0, 0.0]), 0.15)
>  Cylinder(Line([-2.225, -0.7, 0.25], [-1.0, 0.0, 0.0]), 0.15)
> ```
> 
> Rather than the trapezoidal unit cell, we can also specify the boundary as the Wigner-Seitz
> unit cell, using Brillouin.jl's `wignerseitz`:
> ```jl
> julia> using Brillouin
> 
> julia> uc = wignerseitz(Rs);
> 
> julia> uc_cs′ = symmetrize(ops, cs, uc)
>  Cylinder(Line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), 0.15)
>  Cylinder(Line([-0.225, 0.7, 0.25], [-1.0, 0.0, 0.0]), 0.15)
>  Cylinder(Line([-1.7750000000000001, -0.7, -0.25], [-1.0, 0.0, 0.0]), 0.15)
>  Cylinder(Line([-1.7750000000000001, 0.7, -0.25], [-1.0, 0.0, 0.0]), 0.15)
>  Cylinder(Line([-2.225, -0.7, 0.25], [-1.0, 0.0, 0.0]), 0.15)
>  Cylinder(Line([1.55, 0.0, 0.5], [1.0, 0.0, 0.0]), 0.15)
>  Cylinder(Line([-1.55, 0.0, -0.5], [-1.0, 0.0, 0.0]), 0.15)
> ```
> 
> In this case, more cylinders are needed: in general, the number of cylinders may differ for
> trapezoidal and Wigner-Seitz unit cells (but the cylinders' enclosed volume is invariant).
> 
> `Sphere`s work the same way; e.g., for a sphere at a general position of space group 14:
> ```jl
> julia> ss′ = symmetrize(ops, Sphere([0.1, 0.2, 0.05], 0.04), Rs)
>  Sphere([0.1, 0.2, 0.05], 0.04)
>  Sphere([-0.325, -0.5, 0.2], 0.04)
>  Sphere([-0.1, -0.2, -0.05], 0.04)
>  Sphere([0.325, 0.5, -0.2], 0.04)
> ```
> 
> `Cylinder`s and `Sphere`s can also be combined in a single call, by collecting them in a
> `Vector{Primitive}`, e.g., to decorate a woodpile structure with symmetry-related spheres:
> ```jl
> julia> symmetrize(ops, Primitive[cs, Sphere([0.1, 0.2, 0.05], 0.04)], Rs)
> ```
> 
> ## Visualization
> The "symmetrized" primitives can be visualized in the unit cell associated with `boundary`
> using Makie.jl. For example, the following illustrates a symmetric woodpile structure in
> space group 14 in its trapezoidal and Wigner Seitz unit cells, respectively:
> 
> ```jl
> julia> using GLMakie
> 
> julia> plot(cs′, Rs) # trapezoidal unit cell
> 
> julia> plot(uc_cs′, uc) # Wigner-Seitz unit cell
> 
> julia> plot(ss′, Rs) # spheres plot the same way
> ```
> 
> The number of samples used to resolve the primitives' isosurfaces can be controlled using
> the keyword `samples` in the `plot` function. Similarly, the isosurfaces can be plotted
> as `:merged` or `:individual` (default: `:merged`) via the `style` keyword argument.