"""
Post-processing of mesh files
Pre-processing of PD programs
Handle with 2D linear meshes
Processing of standard 3D linear mesh files
"""
mutable struct Post2D <: AbstractPDGeometry
    bc_nodes::Vector{Int}
    pos::Matrix{Float64}
    area::Vector{Float64}
    in_bc_nodes::Vector{Int}
    bc_edge::Dict{Int, Vector{Vector{Vector{Float64}}}}
    mesh_nodes::Dict{Int, Vector{Float64}}
    pd_id_to_surface::Dict{Int, Int}    
    surface_elements::Dict{Int, Vector{Int}}  # Added to store surface elements
    element_type::Dict{Int, String}
    bc_counter::BoundaryCounter
    
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
            nodes, out_boundary_elements, surface_elements, element_type = get_2d_boundary_elements(full_filepath, outside_elsets)
            in_boundary_elements = Dict{Int, Vector{Int}}()
        else
            nodes, out_boundary_elements, in_boundary_elements, surface_elements, element_type = 
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
        pd_id_to_surface = Dict(pd_id => element_id for (element_id, pd_id) in mesh_node_ids)      
        
        # Convert area to vector for consistency
        area_vec = vec(area)
        
        bc_counter = BoundaryCounter()
        
        for (mesh_id, elem) in surface_elements
            pd_id = mesh_node_ids[mesh_id]
            cell_type = element_type[mesh_id]
            add_element!(bc_counter, elem, cell_type, pd_id)
        end
        
        new(bc_nodes, pos, area_vec, in_bc_nodes, bc_edge, nodes, pd_id_to_surface, surface_elements, element_type, bc_counter)
    end
end

# Helper function to access the main components
function Base.show(io::IO, geo::Post2D)
    println(io, "Post2D Structure:")
    println(io, "  - Number of PD nodes: ", size(geo.pos, 2))
    println(io, "  - Number of BC nodes: ", length(geo.bc_nodes))
    println(io, "  - Number of inner BC nodes: ", length(geo.in_bc_nodes))
    println(io, "  - Position matrix size: ", size(geo.pos))
    println(io, "  - Area vector length: ", length(geo.area))
end

# return node Dict, index and position 
#boundary_elements Dict, index and node index
#surface_elements Dict, index and node index


function get_2d_boundary_elements(filepath::String, out_target_elsets::Vector{String})
    out_boundary_elements = Dict{Int, Vector{Int}}()
    surface_elements = Dict{Int, Vector{Int}}()
    element_type = Dict{Int, String}()
    nodes = Dict{Int, Vector{Float64}}()

    current_section = ""
    current_elset = ""
    current_element_type = ""

    if !isfile(filepath)
        error("File does not exist: $filepath")
    end

    for line in eachline(filepath)
        line = strip(line)
        if isempty(line) || startswith(line, "**")
            continue
        end

        # Section headers
        if startswith(line, "*NODE")
            current_section = "Node"
            continue
        elseif startswith(line, "*ELEMENT")
            current_section = "Element"
            # Element type
            m_type = match(r"TYPE=([\w\d]+)"i, line)
            current_element_type = m_type !== nothing ? m_type.captures[1] : ""
            # ELSET
            m_set = match(r"ELSET=([\w\d]+)"i, line)
            current_elset = m_set !== nothing ? m_set.captures[1] : ""
            continue
        elseif startswith(line, "*")
            current_section = ""
            current_elset = ""
            current_element_type = ""
            continue
        end

        # Parse nodes
        if current_section == "Node"
            parts = split(line, ",")
            try
                node_id = parse(Int, strip(parts[1]))
                coord = [parse(Float64, strip(x)) for x in parts[2:end]]
                nodes[node_id] = coord
            catch
                @warn "Invalid node line: $line"
            end

        # Parse elements
        elseif current_section == "Element"
            parts = split(line, ",")
            try
                elem_id = parse(Int, strip(parts[1]))
                elem_nodes = [parse(Int, strip(x)) for x in parts[2:end]]

                # Boundary elements
                if current_elset in out_target_elsets
                    out_boundary_elements[elem_id] = elem_nodes

                # Only surface elements
                elseif occursin("Surface", current_elset)
                    surface_elements[elem_id] = elem_nodes
                    element_type[elem_id] = current_element_type
                end
            catch
                @warn "Invalid element line: $line"
            end
        end
    end

    return nodes, out_boundary_elements, surface_elements, element_type
end

function get_2d_boundary_elements(filepath::String, out_target_elsets::Vector{String}, in_target_elsets::Vector{String})
    out_boundary_elements = Dict{Int, Vector{Int}}()
    in_boundary_elements = Dict{Int, Vector{Int}}()
    surface_elements = Dict{Int, Vector{Int}}()
    element_type = Dict{Int, String}()
    nodes = Dict{Int, Vector{Float64}}()

    current_section = ""
    current_elset = ""
    current_element_type = ""

    if !isfile(filepath)
        error("File does not exist: $filepath")
    end

    for line in eachline(filepath)
        line = strip(line)
        if isempty(line) || startswith(line, "**")
            continue
        end

        # Section headers
        if startswith(line, "*NODE")
            current_section = "Node"
            continue
        elseif startswith(line, "*ELEMENT")
            current_section = "Element"
            # Element type
            m_type = match(r"TYPE=([\w\d]+)"i, line)
            current_element_type = m_type !== nothing ? m_type.captures[1] : ""
            # ELSET
            m_set = match(r"ELSET=([\w\d]+)"i, line)
            current_elset = m_set !== nothing ? m_set.captures[1] : ""
            continue
        elseif startswith(line, "*")
            current_section = ""
            current_elset = ""
            current_element_type = ""
            continue
        end

        # Parse nodes
        if current_section == "Node"
            parts = split(line, ",")
            try
                node_id = parse(Int, strip(parts[1]))
                coord = [parse(Float64, strip(x)) for x in parts[2:end]]
                nodes[node_id] = coord
            catch
                @warn "Invalid node line: $line"
            end

        # Parse elements
        elseif current_section == "Element"
            parts = split(line, ",")
            if length(parts) < 2
                @warn "Invalid element line: $line"
                continue
            end
            try
                elem_id = parse(Int, strip(parts[1]))
                elem_nodes = [parse(Int, strip(x)) for x in parts[2:end]]

                # Out boundary
                if current_elset in out_target_elsets
                    out_boundary_elements[elem_id] = elem_nodes
                # In boundary
                elseif current_elset in in_target_elsets
                    in_boundary_elements[elem_id] = elem_nodes
                # Surface elements
                elseif occursin("Surface", current_elset)
                    surface_elements[elem_id] = elem_nodes
                    element_type[elem_id] = current_element_type
                end
            catch
                @warn "Failed parsing element line: $line"
            end
        end
    end

    return nodes, out_boundary_elements, in_boundary_elements, surface_elements, element_type
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


