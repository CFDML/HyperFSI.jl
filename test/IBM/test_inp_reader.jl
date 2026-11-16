@testitem "read_inp Abaqus part/instance" begin
    using HyperFSI: read_inp, Post2D, pd_sets, section_materials
    mesh = read_inp(joinpath(@__DIR__, "..", "test_data", "abaqus_part_2mat.inp"); verbose=false)
    @test length(mesh.nodes) == 6
    @test length(mesh.elements) == 2
    @test mesh.nodes[1] ≈ [1.0, 0.0, 0.0]          # instance translation applied
    @test mesh.elsets["Steel"] == [1]
    @test mesh.elsets["Copper"] == [2]              # GENERATE
    @test mesh.elsets["All"] == [1, 2]              # set made of other sets
    @test mesh.elsets["Part-1-1.Steel"] == [1]      # qualified name kept
    @test mesh.nsets["Fix"] == [3, 6]               # assembly set with instance=
    @test section_materials(mesh)["Steel"] == "STEEL"

    geo = Post2D(mesh; vtk=nothing)                 # boundary = all free edges
    @test size(geo.pos, 2) == 2
    @test geo.area ≈ [1.0, 1.0]
    @test geo.pos[:, 1] ≈ [1.5, 0.5, 0.0]
    @test sort(geo.bc_nodes) == [1, 2]
    @test sum(length, values(geo.bc_edge)) == 6

    sets = pd_sets(mesh, geo)
    @test sets[:Steel] == [1]
    @test sets[:Copper] == [2]
    @test sets[:Wall] == [1, 2]                     # *SURFACE elset, S1
    @test sets[:Top] == [1, 2]                      # assembly surface
    @test sets[:LeftEdge] == [1]                    # NSET -> elements with a face on it
    @test sets[:Fix] == [2]

    # outer boundary selected by group names (bottom + right + top + left)
    geo2 = Post2D(mesh; boundary=["Wall", "Fix", "Top", "LeftEdge"], inner_boundary=["Fix"], vtk=nothing)
    @test sum(length, values(geo2.bc_edge)) == 6
    @test geo2.in_bc_nodes == [2]
    # an open boundary is rejected unless check=false
    @test_throws ErrorException Post2D(mesh; boundary=["Wall"], vtk=nothing)
    @test sum(length, values(Post2D(mesh; boundary=["Wall"], check=false, vtk=nothing).bc_edge)) == 2
    # domain restricted to one material
    geo3 = Post2D(mesh; domain=["Copper"], vtk=nothing)
    @test size(geo3.pos, 2) == 1
    @test pd_sets(mesh, geo3; domain=["Copper"])[:Copper] == [1]
    # helpful error for unknown group names
    err = try Post2D(mesh; boundary=["Nope"], vtk=nothing) catch e e end
    @test err isa ErrorException && occursin("Available groups", err.msg)
end

@testitem "Post2D outward orientation with a hole" begin
    using HyperFSI: read_inp, Post2D
    mesh = read_inp(joinpath(@__DIR__, "..", "test_data", "Ring.inp"); verbose=false)
    geo = Post2D(mesh; vtk=nothing)                 # outer circle + hole = 2 loops
    # Green's theorem: with outward normals the boundary integral equals the area
    a = sum(0.5 * (s[1][1] * s[2][2] - s[2][1] * s[1][2]) for segs in values(geo.bc_edge) for s in segs)
    @test a ≈ sum(geo.area) rtol=1e-8
    # every segment normal [dy,-dx] points away from its own PD point
    outward(pid, s) = (d = s[2] .- s[1]; m = 0.5 .* (s[1] .+ s[2]); c = geo.pos[:, pid];
                       d[2] * (m[1] - c[1]) - d[1] * (m[2] - c[2]) > 0)
    @test all(outward(pid, s) for (pid, segs) in geo.bc_edge for s in segs)
end

@testitem "Post3D non-convex body and quadratic elements" begin
    using HyperFSI: read_inp, Post3D, pd_sets
    using HyperFSI.LinearAlgebra: cross, dot   # LinearAlgebra is not a direct test dependency
    dir = mktempdir()
    # L-shaped body of three unit hexahedra (non-convex), lower-case keywords,
    # element connectivity split over two lines
    file = joinpath(dir, "lshape.inp")
    cubes = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
    nid = Dict{NTuple{3,Int},Int}()
    io = IOBuffer()
    println(io, "*heading\n** comment line\n*node")
    for (ox, oy, oz) in cubes, z in 0:1, y in 0:1, x in 0:1
        p = (ox + x, oy + y, oz + z)
        haskey(nid, p) && continue
        nid[p] = length(nid) + 1
        println(io, nid[p], ", ", p[1], ", ", p[2], ", ", p[3])
    end
    println(io, "*element, type=C3D8, elset=Body")
    for (k, (ox, oy, oz)) in enumerate(cubes)
        c = [nid[(ox, oy, oz)], nid[(ox+1, oy, oz)], nid[(ox+1, oy+1, oz)], nid[(ox, oy+1, oz)],
             nid[(ox, oy, oz+1)], nid[(ox+1, oy, oz+1)], nid[(ox+1, oy+1, oz+1)], nid[(ox, oy+1, oz+1)]]
        println(io, k, ", ", join(c[1:5], ", "), ",")
        println(io, join(c[6:8], ", "))
    end
    println(io, "*nset, nset=Bottom")
    println(io, join([v for (p, v) in nid if p[3] == 0], ", "))
    write(file, String(take!(io)))

    mesh = read_inp(file; verbose=false)
    geo = Post3D(mesh; vtk=nothing)
    @test size(geo.pos, 2) == 3
    @test geo.vol ≈ ones(3)
    @test sum(length, values(geo.bc_surface)) == 14
    # divergence theorem: (1/3) Σ x_f · A_f = V only if all normals point outward
    flux(f) = dot(sum(f) ./ length(f), 0.5 .* sum(cross(f[i], f[mod1(i + 1, length(f))]) for i in eachindex(f))) / 3
    V = sum(flux(f) for faces in values(geo.bc_surface) for f in faces)
    @test V ≈ 3.0
    @test sort(pd_sets(mesh, geo)[:Bottom]) == [1, 2, 3]
    geo_b = Post3D(mesh; boundary=["Bottom"], check=false, vtk=nothing)
    @test sum(length, values(geo_b.bc_surface)) == 3

    # one quadratic hexahedron (C3D20) written over two lines
    file2 = joinpath(dir, "hex20.inp")
    corners = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
    mids = [(0.5,0,0),(1,0.5,0),(0.5,1,0),(0,0.5,0),(0.5,0,1),(1,0.5,1),(0.5,1,1),(0,0.5,1),
            (0,0,0.5),(1,0,0.5),(1,1,0.5),(0,1,0.5)]
    open(file2, "w") do f
        println(f, "*NODE")
        for (i, p) in enumerate(vcat(corners, mids))
            println(f, i, ", ", p[1], ", ", p[2], ", ", p[3])
        end
        println(f, "*ELEMENT, TYPE=C3D20R, ELSET=EALL")
        println(f, "1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15")
        println(f, "16, 17, 18, 19, 20")
        println(f, "*SURFACE, NAME=TOPSURF\nEALL, S2")
    end
    geo2 = Post3D(read_inp(file2; verbose=false); vtk=nothing)
    @test geo2.vol ≈ [1.0]
    @test geo2.pos[:, 1] ≈ [0.5, 0.5, 0.5]
    @test length(geo2.bc_surface[1]) == 6
    top = Post3D(file2; boundary="TOPSURF", check=false, vtk=nothing)
    @test all(p -> p[3] ≈ 1.0, top.bc_surface[1][1])
end
