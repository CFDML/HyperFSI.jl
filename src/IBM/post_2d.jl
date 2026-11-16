"""
    Post2D(mesh::InpMesh; boundary=:all, inner_boundary=String[], domain=nothing,
           check=true, vtk="boundary_corrected.vtk")
    Post2D(filepath; kwargs...)
    Post2D(filepath, outside_elsets, inside_elsets=String[])      # legacy form

2D pre-processing: every 2D element of the mesh becomes one PD point (centroid, area),
and the boundary used by the IBM coupling is selected by **physical group names**.

Keyword arguments
- `boundary`: outer (fluid-wetted) boundary. Names of physical groups (ELSET of line
  elements, ELSET of 2D elements -> their free edges, SURFACE, NSET), or `:all` for
  every free edge of the mesh (default).
- `inner_boundary`: groups whose adjacent PD points are returned in `in_bc_nodes`.
- `domain`: element sets that form the solid (default: all 2D elements).
- `check`: verify that `boundary` forms closed, consistently oriented loops.
- `vtk`: file to write the oriented boundary to (`nothing` to skip).

Use [`pd_sets`](@ref) to get the PD points of each group (materials, BC regions).

```julia
mesh = read_inp("probe.inp"); println(mesh)      # list the groups of the file
geo  = Post2D(mesh; boundary=["Wall"])
sets = pd_sets(mesh, geo)                          # Dict(:Steel => [...], :Wall => [...])
```
"""
struct Post2D
    bc_nodes::Vector{Int}
    pos::Matrix{Float64}
    area::Vector{Float64}
    in_bc_nodes::Vector{Int}
    bc_edge::Dict{Int, Vector{Vector{Vector{Float64}}}}

    Post2D(bc_nodes, pos, area, in_bc_nodes, bc_edge) = new(bc_nodes, pos, area, in_bc_nodes, bc_edge)
end

function Post2D(mesh::InpMesh; boundary=:all, inner_boundary=String[], domain=nothing,
                check::Bool=true, vtk::Union{Nothing,AbstractString}="boundary_corrected.vtk")
    println("=== Post2D Processing ===")
    dom = domain_elements(mesh, 2; domain=domain)
    pos, area = _pd_points(mesh, dom, 2)
    pd_of = Dict(e => i for (i, e) in enumerate(dom))
    println("PD points: $(length(dom)) (2D elements", domain === nothing ? "" : " of $(domain)", ")")
    println("Outer boundary: ", boundary === :all ? "all free edges" : join(string.(boundary isa AbstractVector ? boundary : [boundary]), ", "))

    faces = _boundary_faces(mesh, dom, boundary, 2)
    isempty(faces) && error("The outer boundary $(boundary) is empty.")
    bc_edge = Dict{Int, Vector{Vector{Vector{Float64}}}}()
    face_nodes = Tuple{Int,Vector{Int}}[]
    for (e, f) in faces
        pid = pd_of[e]
        push!(get!(bc_edge, pid, Vector{Vector{Float64}}[]), [copy(mesh.nodes[f[1]]), copy(mesh.nodes[f[2]])])
        push!(face_nodes, (pid, copy(f)))
    end
    bc_nodes = sort!(collect(keys(bc_edge)))

    in_bc_nodes = Int[]
    if !(inner_boundary === nothing || isempty(inner_boundary))
        in_faces = _boundary_faces(mesh, dom, inner_boundary, 2)
        in_bc_nodes = unique!(sort!([pd_of[e] for (e, _) in in_faces]))
    end

    orient_boundary_2d!(bc_edge, pos)
    check && check_boundary_closed(face_nodes, 2)
    vtk === nothing || export_boundary_vtk(bc_edge, vtk)

    geo = Post2D(bc_nodes, pos, area, in_bc_nodes, bc_edge)
    return geo
end

Post2D(filepath::AbstractString; kwargs...) = Post2D(read_inp(filepath); kwargs...)

# Legacy interface: Post2D("mesh.inp", ["Line1", "Line2"], ["Line5"])
function Post2D(filepath::String, outside_elsets::Vector{String}, inside_elsets::Vector{String}=String[])
    isempty(filepath) && error("File path cannot be empty!")
    isempty(outside_elsets) && error("outside_elsets cannot be empty!")
    return Post2D(read_inp(filepath); boundary=outside_elsets, inner_boundary=inside_elsets)
end

function Base.show(io::IO, geo::Post2D)
    println(io, "Post2D Structure:")
    println(io, "  - Number of PD nodes: ", size(geo.pos, 2))
    println(io, "  - Number of GC nodes: ", length(geo.bc_nodes))
    println(io, "  - Number of inner BC nodes: ", length(geo.in_bc_nodes))
    println(io, "  - Position matrix size: ", size(geo.pos))
    println(io, "  - Area vector length: ", length(geo.area))
end

"""
    orient_boundary_2d!(bc_edge, pos)

Orient every boundary segment so that its outward normal `[dy, -dx]` points away from
the centroid of the element (PD point) it belongs to. This is exact for any shape,
including non-convex bodies and bodies with holes (outer loops become counter-clockwise,
hole loops clockwise).
"""
function orient_boundary_2d!(bc_edge::Dict{Int, Vector{Vector{Vector{Float64}}}}, pos::AbstractMatrix)
    flipped = 0
    for (pid, segs) in bc_edge
        c = pos[:, pid]
        for (k, seg) in enumerate(segs)
            d = seg[2] .- seg[1]
            mid = 0.5 .* (seg[1] .+ seg[2])
            if d[2] * (mid[1] - c[1]) - d[1] * (mid[2] - c[2]) < 0
                segs[k] = [seg[2], seg[1]]
                flipped += 1
            end
        end
    end
    println("Boundary orientation: $flipped of $(sum(length, values(bc_edge); init=0)) segments reversed (outward normals).")
    return flipped
end

"""
    check_boundary_closed(faces, dim)

`faces` are `(pd_id, corner node ids)` pairs. In 2D every node must belong to exactly
two boundary segments (closed loops); in 3D every edge must be shared by exactly two
boundary faces (closed surface). Throws a descriptive error otherwise.
"""
function check_boundary_closed(faces::Vector{Tuple{Int,Vector{Int}}}, dim::Int)
    println("\n=== Checking boundary watertightness ===")
    if dim == 2
        deg = Dict{Int,Int}()
        for (_, f) in faces, n in f
            deg[n] = get(deg, n, 0) + 1
        end
        bad = [n for (n, d) in deg if d != 2]
        if !isempty(bad)
            error("GEOMETRY ERROR: the boundary is not closed: $(length(bad)) node(s) belong to " *
                  "$(join(unique(deg[n] for n in bad), "/")) boundary segment(s) instead of 2 " *
                  "(e.g. node $(first(bad))). Check the boundary groups (missing line, or an internal line selected).")
        end
        # count loops
        adj = Dict{Int,Vector{Int}}()
        for (_, f) in faces
            push!(get!(adj, f[1], Int[]), f[2]); push!(get!(adj, f[2], Int[]), f[1])
        end
        seen = Set{Int}(); loops = 0
        for s in keys(adj)
            s in seen && continue
            loops += 1
            stack = [s]
            while !isempty(stack)
                v = pop!(stack); v in seen && continue; push!(seen, v); append!(stack, adj[v])
            end
        end
        println("✓ Closed boundary: $(length(faces)) segments in $loops loop(s)")
    else
        edges = Dict{Tuple{Int,Int},Int}()
        for (_, f) in faces
            n = length(f)
            for i in 1:n
                a, b = f[i], f[mod1(i + 1, n)]
                key = a < b ? (a, b) : (b, a)
                edges[key] = get(edges, key, 0) + 1
            end
        end
        bad = [k for (k, c) in edges if c != 2]
        if !isempty(bad)
            error("GEOMETRY ERROR: the boundary surface is not closed: $(length(bad)) edge(s) are not shared by " *
                  "exactly 2 boundary faces (e.g. edge $(first(bad))). Check the boundary groups.")
        end
        println("✓ Closed boundary surface: $(length(faces)) faces")
    end
    return true
end

function calculate_normal_vector(node1::Vector{Float64}, node2::Vector{Float64})
    direction = node2 .- node1
    normal = [-direction[2], direction[1], 0.0]
    norm_length = sqrt(sum(normal.^2)) 
    if norm_length == 0
        error("The two points are identical, unable to calculate a valid normal vector.")
    end  
    return normal ./ norm_length
end

function calculate_normal_vector(node1::Vector{Float64}, node2::Vector{Float64}, node3::Vector{Float64})
    v1 = node2 .- node1
    v2 = node3 .- node1
    normal_vector = cross(v1, v2)
    magnitude = norm(normal_vector)
    if magnitude ≈ 0.0
        error("The three points are collinear; a normal vector cannot be defined.")
    end
    unit_normal_vector = normal_vector ./ magnitude
    return unit_normal_vector
end


#Export boundary to VTK format with outward normal vectors.
#Open in ParaView and use 'Glyph' filter to visualize normals.
function export_boundary_vtk(bc_edge::Dict{Int, Vector{Vector{Vector{Float64}}}}, 
                             filename::String="boundary.vtk")
    println("Exporting boundary to VTK format...")
    
    # Collect unique points and segments
    point_map = Dict{Tuple{Float64, Float64, Float64}, Int}()
    points = Vector{Float64}[]
    lines = Tuple{Int, Int}[]
    normals = Vector{Float64}[]
    
    for (pd_id, segments) in bc_edge
        for segment in segments
            length(segment) != 2 && continue
            
            indices = Int[]
            for point in segment
                key = (point[1], point[2], point[3])
                if !haskey(point_map, key)
                    point_map[key] = length(points)
                    push!(points, point)
                end
                push!(indices, point_map[key])
            end
            push!(lines, (indices[1], indices[2]))
            
            # Calculate outward normal (90° CW rotation for CCW boundary)
            dir = segment[2] .- segment[1]
            normal = [dir[2], -dir[1], 0.0]
            len = sqrt(normal[1]^2 + normal[2]^2)
            len > 1e-10 && (normal ./= len)
            push!(normals, normal)
        end
    end
    
    # Write VTK file
    open(filename, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "Boundary edges with normals")
        println(io, "ASCII")
        println(io, "DATASET POLYDATA")
        println(io, "POINTS $(length(points)) float")
        
        for p in points
            println(io, "$(p[1]) $(p[2]) $(p[3])")
        end
        
        println(io, "\nLINES $(length(lines)) $(length(lines) * 3)")
        for (i, j) in lines
            println(io, "2 $i $j")
        end
        
        # Add normal vectors as CELL_DATA
        println(io, "\nCELL_DATA $(length(lines))")
        println(io, "VECTORS normal float")
        for n in normals
            println(io, "$(n[1]) $(n[2]) $(n[3])")
        end
    end
    
    println("✓ VTK file saved to: $filename")
    println("  Open in ParaView -> Filters -> Glyph -> Vectors: 'normal', Glyph Type: 'Arrow'")
end
