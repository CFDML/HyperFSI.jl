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



@testitem "Post2D keyword interface" begin
    using HyperFSI: Post2D, read_inp, pd_sets
    mesh_file = joinpath(@__DIR__, "..", "test_data", "squ.inp")
    geo_all = Post2D(mesh_file; vtk=nothing)                       # all free edges
    geo_old = Post2D(mesh_file, ["Line1", "Line2", "Line3", "Line4"])
    @test geo_all.pos == geo_old.pos
    @test geo_all.bc_nodes == geo_old.bc_nodes
    @test keys(geo_all.bc_edge) == keys(geo_old.bc_edge)
    @test all(Set(geo_all.bc_edge[k]) == Set(geo_old.bc_edge[k]) for k in keys(geo_old.bc_edge))
    mesh = read_inp(mesh_file; verbose=false)
    sets = pd_sets(mesh, geo_all)
    @test sets[:Surface1] == [1, 2, 3, 4]
    @test length(sets[:Line1]) == 2
    @test all(i -> geo_all.pos[2, i] < 0, sets[:Line1])
end
