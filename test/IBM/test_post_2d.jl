@testitem "Post2D" begin
    using  HyperFSI: Post2D
    mesh_file = joinpath(@__DIR__, "..", "test_data", "squ.inp")
    geo = Post2D(mesh_file, ["Line1", "Line2", "Line3", "Line4"])
    pos0 = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]

    @test geo isa Post2D
    @test length(geo.bc_nodes) > 0
    @test sort(geo.bc_nodes) == sort([1, 2, 3, 4])
    @test size(geo.pos, 1) > 0
    @test size(geo.pos, 2) > 0
    @test Main.are_cols_approx_equal(geo.pos, pos0) 
    @test abs(geo.area[1] - 0.25) <= 1e-9
    @test size(geo.pos, 2) == length(geo.area)

    @test_throws ErrorException Post2D("", ["OuterBoundary"])
end


