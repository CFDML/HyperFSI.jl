@testitem "Post3D" begin
    using  HyperFSI: Post3D
    mesh_file = "/Users/shiweihu/Codes/Develops/HyperFSI/test/test_data/vol.inp"
    geo = Post3D(mesh_file, ["Surface1", "Surface2", "Surface3", "Surface4", "Surface5", "Surface6"])
    pos0 = [ -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
             -0.25  -0.25  -0.25  -0.25   0.25   0.25   0.25  0.25]

    @test geo isa Post3D
    @test length(geo.bc_nodes) > 0
    @test sort(geo.bc_nodes) == sort([1, 2, 3, 4, 5, 6, 7, 8])
    @test size(geo.pos, 1) > 0
    @test size(geo.pos, 2) > 0
    @test Main.are_cols_approx_equal(geo.pos, pos0) 
    @test abs(geo.vol[1] - 0.125) <= 1e-9
    @test size(geo.pos, 2) == length(geo.vol)

    @test_throws ErrorException Post3D("", ["OuterBoundary"])
end

