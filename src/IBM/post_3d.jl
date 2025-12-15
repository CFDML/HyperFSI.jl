struct Post3D
    bc_nodes::Vector{Int}
    pos::Matrix{Float64}
    vol::Vector{Float64}
    in_bc_nodes::Vector{Int}
    bc_surface::Dict{Int, Vector{Vector{Vector{Float64}}}}
    
    function Post3D(filepath::String, outside_elsets::Vector{String}, inside_elsets::Vector{String}=String[])
        # Input validation
        isempty(filepath) && error("File path cannot be empty!")
        isempty(outside_elsets) && error("outside_elsets cannot be empty!")
        
        println("=== Post3D Processing ===")
        full_filepath = joinpath(pwd(), filepath) 
        println("Loading file: $full_filepath") 
        println("Processing outer-boundaries: ", join(outside_elsets, ", "))
        
        # Process mesh data
        if isempty(inside_elsets)
            nodes, out_boundary_elements, volume_elements = get_3d_boundary_elements(full_filepath, outside_elsets)
            in_boundary_elements = Dict{Int, Vector{Int}}()
        else
            nodes, out_boundary_elements, in_boundary_elements, volume_elements = 
                get_3d_boundary_elements(full_filepath, outside_elsets, inside_elsets)
        end
        
        node_to_volume = build_node_to_volume_map(volume_elements)
        gc_elements = find_volume_elements_containing_boundary_with_coordinates(
            out_boundary_elements, volume_elements, node_to_volume, nodes)
        
        # Generate PD nodes
        pos, vol, mesh_node_ids = mesh_to_nodes_3d(nodes, volume_elements)
        @assert size(pos, 2) == length(volume_elements) "Mismatch between elements and PD nodes."
        
        # Process boundary nodes
        bc_nodes = Int[]
        bc_surface = Dict{Int, Vector{Vector{Vector{Float64}}}}()
        for element_id in keys(gc_elements)
            pd_id = mesh_node_ids[element_id]
            push!(bc_nodes, pd_id)
            bc_surface[pd_id] = gc_elements[element_id]
        end
        
        # Process inner boundary if exists
        in_bc_nodes = Int[]
        if !isempty(in_boundary_elements)
            in_elements = find_volume_elements_containing_boundary_with_coordinates(
                in_boundary_elements, volume_elements, node_to_volume, nodes)
            for element_id in keys(in_elements)
                pd_id = mesh_node_ids[element_id]
                push!(in_bc_nodes, pd_id)
            end
        end        
        
        # Convert vol to vector for consistency
        vol_vec = vec(vol)

        # Check and fix boundary watertightness
        if !isempty(bc_surface)
            check_and_fix_boundary_watertightness_3d!(bc_surface)
        end
        
        new(bc_nodes, pos, vol_vec, in_bc_nodes, bc_surface)
    end
end

# Helper function to access the main components
function Base.show(io::IO, geo::Post3D)
    println(io, "Post3D Structure:")
    println(io, "  - Number of PD nodes: ", size(geo.pos, 2))
    println(io, "  - Number of GC nodes: ", length(geo.bc_nodes))
    println(io, "  - Number of inner BC nodes: ", length(geo.in_bc_nodes))
    println(io, "  - Position matrix size: ", size(geo.pos))
    println(io, "  - Volume vector length: ", length(geo.vol))
end

# return node Dict, index and position 
#boundary_elements Dict, index and node index
#volume_elements Dict, index and node index

function get_3d_boundary_elements(filepath::String, out_target_elsets::Vector{String}, in_target_elsets::Vector{String})
    out_boundary_elements = Dict{Int, Vector{Int}}()  # Out_boundary surface elements
    in_boundary_elements = Dict{Int, Vector{Int}}()  # In_boundary surface elements    
    volume_elements = Dict{Int, Vector{Int}}()   # Voleme elements
    nodes = Dict{Int, Vector{Float64}}()           # Nodes

    current_set = ""
    current_elset = ""

    println("Reading file: $filepath")
    if !isfile(filepath)
        println("File does not exist!")
        return
    end

    println("Parsing file...")
    for line in eachline(filepath)
        line = strip(line)  # Remove whitespace

        # Parse node section
        if occursin("*NODE", line)
            current_set = "Node"
            continue
        # Parse element section
        elseif occursin("*ELEMENT", line)
            match_result = match(r"\*ELEMENT.*ELSET=(\w+)", line)
            if match_result !== nothing
                current_elset = match_result.captures[1]
                current_set = current_elset
            end
            continue
        # Reset current section for other lines starting with '*'
        elseif startswith(line, "*")
            current_set = ""
            continue
        end

        # Parse nodes
        if current_set == "Node"
            parts = split(line, ",")  # Split by comma
            if length(parts) >= 4  # Check if node format is correct
                node_id = parse(Int, strip(parts[1]))  # Node ID
                coord = [parse(Float64, strip(parts[i])) for i in 2:4]  # Coordinates
                nodes[node_id] = coord
            else
                println("Warning: Incorrect node format -> $line")
            end
        # Parse boundary line or surface elements
        elseif current_set in out_target_elsets || current_set in in_target_elsets || occursin("Volume", current_set)
            parts = split(line, ",")  # Split by comma
            if length(parts) >= 2  # At least element ID and node IDs
                element_id = parse(Int, strip(parts[1]))  # Element ID
                element_nodes = [parse(Int, strip(parts[i])) for i in 2:length(parts)]  # Element nodes

                # Add to boundary elements if in target sets
                if current_set in out_target_elsets
                    out_boundary_elements[element_id] = element_nodes
                elseif current_set in in_target_elsets
                    in_boundary_elements[element_id] = element_nodes                    
                # Add to surface elements if applicable
                elseif occursin("Volume", current_set)
                    volume_elements[element_id] = element_nodes
                end
            else
                println("Warning: Incorrect element format -> $line")
            end
        end
    end

    return nodes, out_boundary_elements, in_boundary_elements, volume_elements
end

function get_3d_boundary_elements(filepath::String, out_target_elsets::Vector{String})
    out_boundary_elements = Dict{Int, Vector{Int}}()  # Out_boundary surface elements
    volume_elements = Dict{Int, Vector{Int}}()   # Voleme elements
    nodes = Dict{Int, Vector{Float64}}()           # Nodes

    current_set = ""
    current_elset = ""

    println("Reading file: $filepath")
    if !isfile(filepath)
        println("File does not exist!")
        return
    end

    println("Parsing file...")
    for line in eachline(filepath)
        line = strip(line)  # Remove whitespace

        # Parse node section
        if occursin("*NODE", line)
            current_set = "Node"
            continue
        # Parse element section
        elseif occursin("*ELEMENT", line)
            match_result = match(r"\*ELEMENT.*ELSET=(\w+)", line)
            if match_result !== nothing
                current_elset = match_result.captures[1]
                current_set = current_elset
            end
            continue
        # Reset current section for other lines starting with '*'
        elseif startswith(line, "*")
            current_set = ""
            continue
        end

        # Parse nodes
        if current_set == "Node"
            parts = split(line, ",")  # Split by comma
            if length(parts) >= 4  # Check if node format is correct
                node_id = parse(Int, strip(parts[1]))  # Node ID
                coord = [parse(Float64, strip(parts[i])) for i in 2:4]  # Coordinates
                nodes[node_id] = coord
            else
                println("Warning: Incorrect node format -> $line")
            end
        # Parse boundary line or surface elements
        elseif current_set in out_target_elsets ||occursin("Volume", current_set)
            parts = split(line, ",")  # Split by comma
            if length(parts) >= 2  # At least element ID and node IDs
                element_id = parse(Int, strip(parts[1]))  # Element ID
                element_nodes = [parse(Int, strip(parts[i])) for i in 2:length(parts)]  # Element nodes

                # Add to boundary elements if in target sets
                if current_set in out_target_elsets
                    out_boundary_elements[element_id] = element_nodes               
                # Add to surface elements if applicable
                elseif occursin("Volume", current_set)
                    volume_elements[element_id] = element_nodes
                end
            else
                println("Warning: Incorrect element format -> $line")
            end
        end
    end

    return nodes, out_boundary_elements, volume_elements
end

function build_node_to_volume_map(volume_elements::Dict{Int, Vector{Int}})
    node_to_volume = Dict{Int, Vector{Int}}()
    
    for (volume_id, nodes) in volume_elements
        for node in nodes
            if haskey(node_to_volume, node)
                push!(node_to_volume[node], volume_id)
            else
                node_to_volume[node] = [volume_id]
            end
        end
    end
    
    return node_to_volume
end

#return Dict for mesh_id to surface_node_position
function find_volume_elements_containing_boundary_with_coordinates(
    boundary_elements::Dict{Int, Vector{Int}}, 
    volume_elements::Dict{Int, Vector{Int}}, 
    node_to_volume::Dict{Int, Vector{Int}}, 
    nodes::Dict{Int, Vector{Float64}})

    containing_elements = Dict{Int, Vector{Vector{Vector{Float64}}}}()

    for (boundary_id, boundary_nodes) in boundary_elements
        candidate_volumes = Set{Int}()
        
        # Find candidate volumes for the current boundary
        for node in boundary_nodes
            if haskey(node_to_volume, node)
                union!(candidate_volumes, node_to_volume[node])
            end
        end
        
        # Check if boundary_nodes are contained in each candidate volume
        for volume_id in candidate_volumes
            if all(node in volume_elements[volume_id] for node in boundary_nodes)
                # Get the coordinates of the boundary nodes
                boundary_coordinates = [nodes[node] for node in boundary_nodes]
                
                # Add the boundary coordinates to the volume in containing_elements
                if haskey(containing_elements, volume_id)
                    push!(containing_elements[volume_id], boundary_coordinates)
                else
                    containing_elements[volume_id] = [boundary_coordinates]
                end
            end
        end
    end

    return containing_elements
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

#3d mesh to pd_node, return pos, vol, and a Dict for mesh_id to pd_id
function mesh_to_nodes_3d(nodes::Dict{Int, Vector{Float64}}, volume_elements::Dict{Int, Vector{Int}})
    N = length(volume_elements)

    pos = zeros(3,N)
    vol = zeros(1,N)
    mesh_node_ids = Dict{Int64, Int64}()

    for (i, (element_id, node_ids)) in enumerate(volume_elements) 
        nodes_coords = [nodes[node_id] for node_id in node_ids]
        x_c = 0.0
        y_c = 0.0
        z_c = 0.0
        area = 0.0
        n = length(nodes_coords)
        for j in 1:n
            x_c += nodes_coords[j][1]
            y_c += nodes_coords[j][2]
            z_c += nodes_coords[j][3]
        end
        pos[:, i] = [x_c/n, y_c/n, z_c/n] 

        if n == 4
            vol[i] = tetvol(nodes_coords[1], nodes_coords[2], nodes_coords[3], nodes_coords[4])
        elseif n == 8
            v1 = tetvol(nodes_coords[1], nodes_coords[2], nodes_coords[4], nodes_coords[5])
            v2 = tetvol(nodes_coords[2], nodes_coords[3], nodes_coords[4], nodes_coords[7])
            v3 = tetvol(nodes_coords[2], nodes_coords[5], nodes_coords[6], nodes_coords[7])
            v4 = tetvol(nodes_coords[4], nodes_coords[5], nodes_coords[7], nodes_coords[8])
            v5 = tetvol(nodes_coords[2], nodes_coords[4], nodes_coords[5], nodes_coords[7])
            vol[i] = v1 + v2 + v3 + v4 + v5
        elseif n == 6
            v1 = tetvol(nodes_coords[1], nodes_coords[2], nodes_coords[3], nodes_coords[6])
            v2 = tetvol(nodes_coords[1], nodes_coords[2], nodes_coords[6], nodes_coords[5])
            v3 = tetvol(nodes_coords[1], nodes_coords[5], nodes_coords[6], nodes_coords[4])
            vol[i] = v1 + v2 + v3
        else
            throw("Only for tetrahedral, wedge or hexahedral elements, currently")
        end

        mesh_node_ids[element_id] = i
    end

    return pos, vol, mesh_node_ids
end

function tetvol(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, d::Vector{Float64})
    v1 = b - a
    v2 = c - a
    v3 = d - a
    vol = 1 / 6 * abs(dot(v1, cross(v2, v3)))
    return vol
end

#Check and fix boundary watertightness and orientation for 3D surfaces.
#Modifies bc_surface dict values in-place while preserving keys.
function check_and_fix_boundary_watertightness_3d!(bc_surface::Dict{Int, Vector{Vector{Vector{Float64}}}})
    println("\n=== Checking 3D Boundary Watertightness ===")
    
    # Extract all boundary faces
    point_coords = Dict{Tuple{Float64, Float64, Float64}, Int}()
    point_id_to_coords = Dict{Int, Vector{Float64}}()
    point_counter = 0
    
    all_faces = []  # [(pd_id, face_idx, face_vertices)]
    
    for (pd_id, faces) in bc_surface
        for (face_idx, face) in enumerate(faces)
            if length(face) < 3
                @warn "Face $face_idx in pd_id $pd_id has < 3 vertices. Skipping."
                continue
            end
            
            push!(all_faces, (pd_id, face_idx, face))
            
            # Register points
            for point in face
                coord_tuple = (round(point[1], digits=8), round(point[2], digits=8), round(point[3], digits=8))
                if !haskey(point_coords, coord_tuple)
                    point_counter += 1
                    point_coords[coord_tuple] = point_counter
                    point_id_to_coords[point_counter] = point
                end
            end
        end
    end
    
    println("Info: Found $(point_counter) unique points, $(length(all_faces)) faces")
    
    if point_counter < 4
        error("ERROR: Too few points ($(point_counter)) to form a valid 3D boundary!")
    end
    
    # Build edge-to-faces mapping for manifold check
    edge_to_faces = Dict{Tuple{Int, Int}, Vector{Int}}()
    
    for (face_id, (pd_id, face_idx, face)) in enumerate(all_faces)
        n = length(face)
        for i in 1:n
            p1_coord = (round(face[i][1], digits=8), round(face[i][2], digits=8), round(face[i][3], digits=8))
            p2_coord = (round(face[i % n + 1][1], digits=8), round(face[i % n + 1][2], digits=8), round(face[i % n + 1][3], digits=8))
            
            u = point_coords[p1_coord]
            v = point_coords[p2_coord]
            
            # Use sorted edge (undirected)
            edge = u < v ? (u, v) : (v, u)
            
            if !haskey(edge_to_faces, edge)
                edge_to_faces[edge] = Int[]
            end
            push!(edge_to_faces[edge], face_id)
        end
    end
    
    # Check manifold: each edge must be shared by exactly 2 faces
    non_manifold_edges = 0
    for (edge, faces) in edge_to_faces
        if length(faces) != 2
            non_manifold_edges += 1
            if non_manifold_edges <= 5
                p1 = point_id_to_coords[edge[1]]
                p2 = point_id_to_coords[edge[2]]
                @warn "Non-manifold edge: $(edge[1])-$(edge[2]) shared by $(length(faces)) faces"
            end
        end
    end
    
    if non_manifold_edges > 0
        error("GEOMETRY ERROR: Found $non_manifold_edges non-manifold edges! " *
              "Boundary must be a closed manifold surface.")
    end
    
    println("✓ Topology valid: closed manifold surface (each edge shared by 2 faces)")
    
    # Calculate geometric center
    center = zeros(3)
    for coord in values(point_id_to_coords)
        center .+= coord
    end
    center ./= point_counter
    println("Geometric center: [$(center[1]), $(center[2]), $(center[3])]")
    
    # Check and fix face orientations
    fixed_count = 0
    
    for (pd_id, face_idx, face) in all_faces
        n = length(face)
        
        # Calculate face center
        face_center = zeros(3)
        for vertex in face
            face_center .+= vertex
        end
        face_center ./= n
        
        # Calculate face normal using cross product
        v1 = face[2] .- face[1]
        v2 = face[3] .- face[1]
        normal = cross(v1, v2)
        len = norm(normal)
        len > 1e-10 && (normal ./= len)
        
        # Vector from geometric center to face center
        outward_dir = face_center .- center
        
        # Check if normal points outward
        if dot(normal, outward_dir) < 0
            # Normal points inward, flip face vertex order
            reverse!(face)
            bc_surface[pd_id][face_idx] = face
            fixed_count += 1
        end
    end
    
    println("Info: Fixed $fixed_count face orientations")
    
    # Generate visualization
    export_boundary_vtk_3d(bc_surface, "boundary_corrected.vtk")
    
    println("=== 3D Boundary Check Complete ===\n")
    
    return fixed_count
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

