using Woodpile
using Test

@testset "Woodpile.jl" begin
    include("t_transforms.jl")
    include("t_cylinder_polygon_intersect.jl")
    include("t_sphere_polygon_intersect.jl")
    include("t_symmetrize.jl")
    include("t_plotting.jl")
end
