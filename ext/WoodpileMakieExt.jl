module WoodpileMakieExt

import Makie: plot
using Makie: Figure, Axis3, Vec, Rect, mesh!, wireframe!, FigureAxisPlot, GeometryBasics
const GB = GeometryBasics
using Meshing: isosurface
using StaticArrays: SVector

using Woodpile: box_boundary, signed_distance, Cylinder, Sphere

function plot(
    cs::AbstractVector{<:Union{Cylinder, Sphere}};
    style::Symbol=:merged,
    boxtol::Real=1e-8,
    limits::NTuple{6, Float64}=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
    samples::Union{NTuple{3, Int}, Int}=50,
    plot_kwargs...
)
    f = Figure()
    ax = Axis3(f[1,1]; aspect=:data, limits=limits)

    # 3D unit cell
    r_o = SVector(limits[1], limits[3], limits[5])
    r_w = SVector(limits[2]-limits[1], limits[4]-limits[3], limits[6]-limits[5])
    p = wireframe!(Rect(r_o, r_w); color=:black, linewidth=0.5)

    # show cylinders / spheres, cropped to 3D unit cell bounding box
    plot_opts = merge((; color=:gray, transparency=false), plot_kwargs)

    xs = range(limits[1], limits[2], length=samples isa NTuple ? samples[1] : samples)
    ys = range(limits[3], limits[4], length=samples isa NTuple ? samples[2] : samples)
    zs = range(limits[5], limits[6], length=samples isa NTuple ? samples[3] : samples)
    if style == :individual
        for c in cs
            V = [(r=SVector(x,y,z); max(signed_distance(r, c), box_boundary(r, limits, boxtol))) for x in xs, y in ys, z in zs]
            vs, _fs = isosurface(V, xs, ys, zs)
            fs = map(f -> GB.TriangleFace(f[1], f[2], f[3]), _fs)
            mesh!(vs, fs; plot_opts...)
        end
    elseif style == :merged
        V = [(r=SVector(x,y,z); max(signed_distance(r, cs), box_boundary(r, limits, boxtol))) for x in xs, y in ys, z in zs]
        vs, _fs = isosurface(V, xs, ys, zs)
        fs = map(f -> GB.TriangleFace(f[1], f[2], f[3]), _fs)
        mesh!(vs, fs; plot_opts...)
    else
        error("unsupported style $style")
    end

    return FigureAxisPlot(f, ax, p)
end

end # module WoodpileMakieExt