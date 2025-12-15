
struct Post2D
    bc_nodes::Vector{Int}
    pos::Matrix{Float64}
    area::Vector{Float64}
    in_bc_nodes::Vector{Int}
    bc_edge::Dict{Int, Vector{Vector{Vector{Float64}}}}
    
    function Post2D(filepath::String, outside_elsets::Vector{String}, inside_elsets::Vector{String}=String[])
        # Input validation
        isempty(filepath) && error("File path cannot be empty!")
        isempty(outside_elsets) && error("outside_elsets cannot be empty!")
        
        println("=== IBM2D Processing ===")
        full_filepath = joinpath(pwd(), filepath) 
        println("Loading file: $full_filepath") 
        println("Processing outer-boundaries: ", join(outside_elsets, ", "))
        
        # Process mesh data
        if isempty(inside_elsets)
            nodes, out_boundary_elements, surface_elements = get_2d_boundary_elements(full_filepath, outside_elsets)
            in_boundary_elements = Dict{Int, Vector{Int}}()
        else
            nodes, out_boundary_elements, in_boundary_elements, surface_elements = 
                get_2d_boundary_elements(full_filepath, outside_elsets, inside_elsets)
        end
        
        node_to_surface = build_node_to_surface_map(surface_elements)
        gc_elements = find_surface_elements_containing_boundary_with_coordinates(
            out_boundary_elements, surface_elements, node_to_surface, nodes)
        
        # Generate PD nodes
        pos, area, mesh_node_ids = mesh_to_nodes_2d(nodes, surface_elements)
        @assert size(pos, 2) == length(surface_elements) "Mismatch between elements and PD nodes."
        
        # Process boundary nodes
        bc_nodes = Int[]
        bc_edge = Dict{Int, Vector{Vector{Vector{Float64}}}}()
        for element_id in keys(gc_elements)
            pd_id = mesh_node_ids[element_id]
            push!(bc_nodes, pd_id)
            bc_edge[pd_id] = deepcopy(gc_elements[element_id])
        end
        
        # Process inner boundary if exists
        in_bc_nodes = Int[]
        if !isempty(in_boundary_elements)
            in_elements = find_surface_elements_containing_boundary_with_coordinates(
                in_boundary_elements, surface_elements, node_to_surface, nodes)
            for element_id in keys(in_elements)
                pd_id = mesh_node_ids[element_id]
                push!(in_bc_nodes, pd_id)
            end
        end        
        
        # Convert area to vector for consistency
        area_vec = vec(area)

        if !isempty(bc_edge)
            check_and_fix_boundary_watertightness!(bc_edge)
        end
        
        new(bc_nodes, pos, area_vec, in_bc_nodes, bc_edge)
    end
end

# Helper function to access the main components
function Base.show(io::IO, geo::Post2D)
    println(io, "Post2D Structure:")
    println(io, "  - Number of PD nodes: ", size(geo.pos, 2))
    println(io, "  - Number of GC nodes: ", length(geo.bc_nodes))
    println(io, "  - Number of inner BC nodes: ", length(geo.in_bc_nodes))
    println(io, "  - Position matrix size: ", size(geo.pos))
    println(io, "  - Area vector length: ", length(geo.area))
end

# return node Dict, index and position 
#boundary_elements Dict, index and node index
#surface_elements Dict, index and node index

function get_2d_boundary_elements(filepath::String, out_target_elsets::Vector{String})
    out_boundary_elements = Dict{Int, Vector{Int}}()  # Boundary line elements
    surface_elements = Dict{Int, Vector{Int}}()   # Surface elements
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
        elseif current_set in out_target_elsets || occursin("Surface", current_set)
            parts = split(line, ",")  # Split by comma
            if length(parts) >= 2  # At least element ID and node IDs
                element_id = parse(Int, strip(parts[1]))  # Element ID
                element_nodes = [parse(Int, strip(parts[i])) for i in 2:length(parts)]  # Element nodes

                # Add to boundary elements if in target sets
                if current_set in out_target_elsets
                    out_boundary_elements[element_id] = element_nodes
                # Add to surface elements if applicable
                elseif occursin("Surface", current_set)
                    surface_elements[element_id] = element_nodes
                end
            else
                println("Warning: Incorrect element format -> $line")
            end
        end
    end

    return nodes, out_boundary_elements, surface_elements
end

function get_2d_boundary_elements(filepath::String, out_target_elsets::Vector{String}, in_target_elsets::Vector{String})
    out_boundary_elements = Dict{Int, Vector{Int}}()  # Out_Boundary line elements
    in_boundary_elements = Dict{Int, Vector{Int}}()  # Inner_Boundary line elements
    surface_elements = Dict{Int, Vector{Int}}()   # Surface elements
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
        elseif current_set in out_target_elsets || current_set in in_target_elsets || occursin("Surface", current_set)
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
                elseif occursin("Surface", current_set)
                    surface_elements[element_id] = element_nodes
                end
            else
                println("Warning: Incorrect element format -> $line")
            end
        end
    end

    return nodes, out_boundary_elements, in_boundary_elements, surface_elements
end

function build_node_to_surface_map(surface_elements::Dict{Int, Vector{Int}})
    node_to_surface = Dict{Int, Vector{Int}}()
    
    for (surface_id, nodes) in surface_elements
        for node in nodes
            if haskey(node_to_surface, node)
                push!(node_to_surface[node], surface_id)
            else
                node_to_surface[node] = [surface_id]
            end
        end
    end
    
    return node_to_surface
end

#return Dict for mesh_id to edge_node_position
function find_surface_elements_containing_boundary_with_coordinates(
    boundary_elements::Dict{Int, Vector{Int}}, 
    surface_elements::Dict{Int, Vector{Int}}, 
    node_to_surface::Dict{Int, Vector{Int}}, 
    nodes::Dict{Int, Vector{Float64}})

    containing_elements = Dict{Int, Vector{Vector{Vector{Float64}}}}()

    for (boundary_id, boundary_nodes) in boundary_elements
        candidate_surfaces = Set{Int}()
        
        # Find candidate surfaces for the current boundary
        for node in boundary_nodes
            if haskey(node_to_surface, node)
                union!(candidate_surfaces, node_to_surface[node])
            end
        end
        
        # Check if boundary_nodes are contained in each candidate surface
        for surface_id in candidate_surfaces
            if all(node in surface_elements[surface_id] for node in boundary_nodes)
                # Get the coordinates of the boundary nodes
                boundary_coordinates = [nodes[node] for node in boundary_nodes]
                
                # Add the boundary coordinates to the surface in containing_elements
                if haskey(containing_elements, surface_id)
                    push!(containing_elements[surface_id], boundary_coordinates)
                else
                    containing_elements[surface_id] = [boundary_coordinates]
                end
            end
        end
    end

    return containing_elements
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

#2d mesh to pd_node, return pos, area, and a Dict for mesh_id to pd_id
function mesh_to_nodes_2d(nodes::Dict{Int, Vector{Float64}}, surface_elements::Dict{Int, Vector{Int}})
    N = length(surface_elements)

    pos = zeros(3,N)
    area = zeros(1,N)
    mesh_node_ids = Dict{Int64, Int64}()

    for (i, (element_id, node_ids)) in enumerate(surface_elements) 
        nodes_coords = [nodes[node_id] for node_id in node_ids]
        x_c = 0.0
        y_c = 0.0
        z_c = 0.0
        area_c = 0.0
        n = length(nodes_coords)
        for j in 1:n
            x_c += nodes_coords[j][1]
            y_c += nodes_coords[j][2]
            z_c += nodes_coords[j][3]
            k = mod(j, n) + 1 
            area_c += nodes_coords[j][1] * nodes_coords[k][2] - nodes_coords[k][1] * nodes_coords[j][2]
        end
        pos[:, i] = [x_c/n, y_c/n, z_c/n] 
        area[i] = 0.5 * abs(area_c) 
        mesh_node_ids[element_id] = i
    end

    return pos, area, mesh_node_ids
end

# Check and fix boundary watertightness and orientation.
# Modifies bc_edge dict values in-place while preserving keys.
function check_and_fix_boundary_watertightness!(bc_edge::Dict{Int, Vector{Vector{Vector{Float64}}}})
    println("\n=== Checking Boundary Watertightness ===")
    
    # Extract all boundary segments and build point mapping
    point_coords = Dict{Tuple{Float64, Float64}, Int}()
    point_id_to_coords = Dict{Int, Vector{Float64}}()
    point_counter = 0
    
    all_segments = []
    
    for (pd_id, segments) in bc_edge
        for (seg_idx, segment) in enumerate(segments)
            length(segment) != 2 && continue
            
            push!(all_segments, (pd_id, seg_idx, segment))
            
            # Register points (using x,y only, ignore z)
            for point in segment
                coord_tuple = (round(point[1], digits=8), round(point[2], digits=8))
                if !haskey(point_coords, coord_tuple)
                    point_counter += 1
                    point_coords[coord_tuple] = point_counter
                    point_id_to_coords[point_counter] = point
                end
            end
        end
    end
    
    println("Info: Found $(point_counter) unique points, $(length(all_segments)) segments")
    
    if point_counter < 3
        error("ERROR: Too few points ($(point_counter)) to form a valid boundary!")
    end
    
    # Build undirected adjacency graph
    adj = Dict{Int, Vector{Int}}()
    for id in 1:point_counter
        adj[id] = Int[]
    end
    
    for (pd_id, seg_idx, segment) in all_segments
        p1_coord = (round(segment[1][1], digits=8), round(segment[1][2], digits=8))
        p2_coord = (round(segment[2][1], digits=8), round(segment[2][2], digits=8))
        
        u = point_coords[p1_coord]
        v = point_coords[p2_coord]
        
        if !(v in adj[u])
            push!(adj[u], v)
        end
        if !(u in adj[v])
            push!(adj[v], u)
        end
    end
    
    # Check topology: all points must have degree 2
    for (id, neighbors) in adj
        if length(neighbors) != 2
            coord = point_id_to_coords[id]
            error("GEOMETRY ERROR: Point $id at ($(coord[1]), $(coord[2])) has degree $(length(neighbors)). " *
                  "Boundary must be a simple closed loop!")
        end
    end
    
    println("✓ Topology valid: all points have degree 2")
    
    # Trace closed loop
    ordered_points = Int[]
    start_node = 1
    curr = start_node
    prev = -1
    
    for _ in 1:point_counter
        push!(ordered_points, curr)
        neighbors = adj[curr]
        next_node = (neighbors[1] == prev) ? neighbors[2] : neighbors[1]
        prev = curr
        curr = next_node
    end
    
    println("✓ Traced closed loop with $(length(ordered_points)) points")
    
    # Calculate signed area using Shoelace formula (2D projection on xy-plane)
    area = 0.0
    n = length(ordered_points)
    for i in 1:n
        p1 = point_id_to_coords[ordered_points[i]]
        p2 = point_id_to_coords[ordered_points[i % n + 1]]
        area += p1[1] * p2[2] - p2[1] * p1[2]
    end
    area /= 2.0
    
    println("Signed area: $area, Absolute area: $(abs(area))")
    
    # Ensure counter-clockwise orientation
    if area < 0
        println("⚠ Reversing from clockwise to counter-clockwise")
        reverse!(ordered_points)
    else
        println("✓ Orientation is counter-clockwise")
    end
    
    # Build correct directed edge set
    correct_edges = Set{Tuple{Int, Int}}()
    for i in 1:n
        u = ordered_points[i]
        v = ordered_points[i % n + 1]
        push!(correct_edges, (u, v))
    end
    
    # Fix segment directions in bc_edge (in-place modification)
    fixed_count = 0
    
    for (pd_id, seg_idx, segment) in all_segments
        p1_coord = (round(segment[1][1], digits=8), round(segment[1][2], digits=8))
        p2_coord = (round(segment[2][1], digits=8), round(segment[2][2], digits=8))
        
        u = point_coords[p1_coord]
        v = point_coords[p2_coord]
        
        if (u, v) in correct_edges
            # Direction correct, no change needed
        elseif (v, u) in correct_edges
            # Flip direction
            bc_edge[pd_id][seg_idx] = [segment[2], segment[1]]
            fixed_count += 1
        else
            @warn "Segment $seg_idx in pd_id $pd_id is not part of the traced loop!"
        end
    end
    
    println("Info: Fixed $fixed_count segment directions")
    
    # Generate visualization
    export_boundary_vtk(bc_edge, "boundary_corrected.vtk")
    
    println("=== Boundary Check Complete ===\n")
    
    return fixed_count
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
