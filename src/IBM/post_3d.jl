"""
    Post3D(mesh::InpMesh; boundary=:all, inner_boundary=String[], domain=nothing,
           check=true, vtk="boundary_corrected.vtk")
    Post3D(filepath; kwargs...)
    Post3D(filepath, outside_elsets, inside_elsets=String[])      # legacy form

3D pre-processing: every 3D element (tet, pyramid, wedge, hex; linear or quadratic)
becomes one PD point (centroid, volume). Boundaries are selected by physical group
names exactly as in [`Post2D`](@ref) (ELSET of surface elements, ELSET of 3D elements ->
free faces, SURFACE `elset, S1`, NSET, or `:all`).
"""
struct Post3D
    bc_nodes::Vector{Int}
    pos::Matrix{Float64}
    vol::Vector{Float64}
    in_bc_nodes::Vector{Int}
    bc_surface::Dict{Int, Vector{Vector{Vector{Float64}}}}

    Post3D(bc_nodes, pos, vol, in_bc_nodes, bc_surface) = new(bc_nodes, pos, vol, in_bc_nodes, bc_surface)
end

function Post3D(mesh::InpMesh; boundary=:all, inner_boundary=String[], domain=nothing,
                check::Bool=true, vtk::Union{Nothing,AbstractString}="boundary_corrected.vtk")
    println("=== Post3D Processing ===")
    dom = domain_elements(mesh, 3; domain=domain)
    pos, vol = _pd_points(mesh, dom, 3)
    pd_of = Dict(e => i for (i, e) in enumerate(dom))
    println("PD points: $(length(dom)) (3D elements", domain === nothing ? "" : " of $(domain)", ")")
    println("Outer boundary: ", boundary === :all ? "all free faces" : join(string.(boundary isa AbstractVector ? boundary : [boundary]), ", "))

    ctx = _face_ctx(mesh, dom)
    faces = _boundary_faces(mesh, dom, boundary, 3; ctx=ctx)
    isempty(faces) && error("The outer boundary $(boundary) is empty.")
    bc_surface = Dict{Int, Vector{Vector{Vector{Float64}}}}()
    face_nodes = Tuple{Int,Vector{Int}}[]
    for (e, f) in faces
        pid = pd_of[e]
        push!(get!(bc_surface, pid, Vector{Vector{Float64}}[]), [copy(mesh.nodes[n]) for n in f])
        push!(face_nodes, (pid, copy(f)))
    end
    bc_nodes = sort!(collect(keys(bc_surface)))

    in_bc_nodes = Int[]
    if !(inner_boundary === nothing || isempty(inner_boundary))
        in_faces = _boundary_faces(mesh, dom, inner_boundary, 3; ctx=ctx)
        in_bc_nodes = unique!(sort!([pd_of[e] for (e, _) in in_faces]))
    end

    orient_boundary_3d!(bc_surface, pos)
    check && check_boundary_closed(face_nodes, 3)
    vtk === nothing || export_boundary_vtk_3d(bc_surface, vtk)

    geo = Post3D(bc_nodes, pos, vol, in_bc_nodes, bc_surface)
    return geo
end

Post3D(filepath::AbstractString; kwargs...) = Post3D(read_inp(filepath); kwargs...)

# Legacy interface: Post3D("mesh.inp", ["Surface1", "Surface2"], ["Surface7"])
function Post3D(filepath::String, outside_elsets::Vector{String}, inside_elsets::Vector{String}=String[])
    isempty(filepath) && error("File path cannot be empty!")
    isempty(outside_elsets) && error("outside_elsets cannot be empty!")
    return Post3D(read_inp(filepath); boundary=outside_elsets, inner_boundary=inside_elsets)
end

function Base.show(io::IO, geo::Post3D)
    println(io, "Post3D Structure:")
    println(io, "  - Number of PD nodes: ", size(geo.pos, 2))
    println(io, "  - Number of GC nodes: ", length(geo.bc_nodes))
    println(io, "  - Number of inner BC nodes: ", length(geo.in_bc_nodes))
    println(io, "  - Position matrix size: ", size(geo.pos))
    println(io, "  - Volume vector length: ", length(geo.vol))
end

"""
    orient_boundary_3d!(bc_surface, pos)

Orient every boundary face so that its normal `(v2-v1)×(v3-v1)` points away from the
centroid of the element (PD point) it belongs to. Exact for non-convex bodies and
bodies with cavities (unlike a test against the global geometric centre).
"""
function orient_boundary_3d!(bc_surface::Dict{Int, Vector{Vector{Vector{Float64}}}}, pos::AbstractMatrix)
    flipped = 0
    for (pid, faces) in bc_surface
        c = pos[:, pid]
        for face in faces
            fc = sum(face) ./ length(face)
            nrm = zeros(3)
            n = length(face)
            for i in 1:n   # Newell normal, robust for warped quads
                a, b = face[i], face[mod1(i + 1, n)]
                nrm .+= cross(a, b)
            end
            if dot(nrm, fc .- c) < 0
                reverse!(face)
                flipped += 1
            end
        end
    end
    println("Boundary orientation: $flipped of $(sum(length, values(bc_surface); init=0)) faces reversed (outward normals).")
    return flipped
end

function calculate_normal_vector_3d(node1::Vector{Float64}, node2::Vector{Float64}, node3::Vector{Float64})
    da = node2 .- node1
    db = node3 .- node1
    normal = [da[2]db[3]-da[3]db[2], -da[1]db[3]+da[3]db[1], da[1]db[2]-da[2]db[1]]
    norm_length = sqrt(sum(normal.^2)) 
    if norm_length == 0
        error("The three points are collinear, unable to calculate a valid normal vector.")
    end  
    return normal ./ norm_length
end


function tetvol(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, d::Vector{Float64})
    v1 = b - a
    v2 = c - a
    v3 = d - a
    vol = 1 / 6 * abs(dot(v1, cross(v2, v3)))
    return vol
end

#Export 3D boundary surface to VTK format with outward normal vectors.
#Open in ParaView and use 'Glyph' filter to visualize normals.
function export_boundary_vtk_3d(bc_surface::Dict{Int, Vector{Vector{Vector{Float64}}}}, 
                             filename::String="boundary.vtk")
    println("Exporting 3D boundary to VTK format...")
    
    # Collect unique points and faces
    point_map = Dict{Tuple{Float64, Float64, Float64}, Int}()
    points = Vector{Float64}[]
    faces = Vector{Int}[]
    normals = Vector{Float64}[]
    
    for (pd_id, face_list) in bc_surface
        for face in face_list
            length(face) < 3 && continue
            
            indices = Int[]
            for vertex in face
                key = (vertex[1], vertex[2], vertex[3])
                if !haskey(point_map, key)
                    point_map[key] = length(points)
                    push!(points, vertex)
                end
                push!(indices, point_map[key])
            end
            push!(faces, indices)
            
            # Calculate face normal
            v1 = face[2] .- face[1]
            v2 = face[3] .- face[1]
            normal = cross(v1, v2)
            len = norm(normal)
            len > 1e-10 && (normal ./= len)
            push!(normals, normal)
        end
    end
    
    # Write VTK file
    open(filename, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "3D Boundary surface with normals")
        println(io, "ASCII")
        println(io, "DATASET POLYDATA")
        println(io, "POINTS $(length(points)) float")
        
        for p in points
            println(io, "$(p[1]) $(p[2]) $(p[3])")
        end
        
        # Count total indices needed
        total_indices = sum(length(f) + 1 for f in faces)
        
        println(io, "\nPOLYGONS $(length(faces)) $total_indices")
        for f in faces
            print(io, "$(length(f))")
            for idx in f
                print(io, " $idx")
            end
            println(io)
        end
        
        # Add normal vectors as CELL_DATA
        println(io, "\nCELL_DATA $(length(faces))")
        println(io, "VECTORS normal float")
        for n in normals
            println(io, "$(n[1]) $(n[2]) $(n[3])")
        end
    end
    
    println("✓ VTK file saved to: $filename")
    println("  Open in ParaView -> Filters -> Glyph -> Vectors: 'normal', Glyph Type: 'Arrow'")
end

