
function update_new_edges!(geo::Post2D, new_ablation_point_idx::Vector{Int})

    new_ablation_elements = Vector{Vector{Int}}()
    for i in new_ablation_point_idx
        mesh_id = geo.pd_id_to_surface[i]
        push!(new_ablation_elements, geo.surface_elements[mesh_id])
    end

    for (i, elem) in enumerate(new_ablation_elements)     
        pd_id = new_ablation_point_idx[i]  
        cell_type = geo.element_type[geo.pd_id_to_surface[pd_id]]
        remove_element!(geo.bc_counter, elem, cell_type, pd_id)
    end

    boundary_keys = get_boundary(geo.bc_counter)

    new_bc_edges = Dict{Int, Vector{Vector{Vector{Float64}}}}()
    
    for (face, pd_id) in boundary_keys
        coords = [geo.mesh_nodes[n] for n in face]
        if haskey(new_bc_edges, pd_id)
            push!(new_bc_edges[pd_id], coords)
        else
            new_bc_edges[pd_id] = [coords]
        end
    end

    geo.bc_edge = new_bc_edges
end

function update_new_edges!(geo::Post3D, new_ablation_point_idx::Vector{Int})

    new_ablation_elements = Vector{Vector{Int}}()
    for i in new_ablation_point_idx
        mesh_id = geo.pd_id_to_volume[i]
        push!(new_ablation_elements, geo.volume_elements[mesh_id])
    end

    for (i, elem) in enumerate(new_ablation_elements)     
        pd_id = new_ablation_point_idx[i]  
        cell_type = geo.element_type[geo.pd_id_to_volume[pd_id]]
        remove_element!(geo.bc_counter, elem, cell_type, pd_id)
    end

    boundary_keys = get_boundary(geo.bc_counter)

    new_bc_surfaces = Dict{Int, Vector{Vector{Vector{Float64}}}}()
    
    for (face, pd_id) in boundary_keys
        coords = [geo.mesh_nodes[n] for n in face]
        if haskey(new_bc_surfaces, pd_id)
            push!(new_bc_surfaces[pd_id], coords)
        else
            new_bc_surfaces[pd_id] = [coords]
        end
    end

    geo.bc_surface = new_bc_surfaces
end

function save_bc_edges_vtk(geo::Post2D, root::String, step::Int, Δt::Float64; pvd_name="bc_edges.pvd")

    mkpath(root)
    new_bc_edges = geo.bc_edge
    vtkfile =  joinpath(root, "bc_edges_" * lpad(string(step), 4, '0') * ".vtk")

    all_points = Vector{Vector{Float64}}()
    polylines = Vector{Vector{Int}}()
    polyline_bc_ids = Vector{Int}()

    point_count = 0
    for (bc_id, edges) in new_bc_edges
        for edge_points in edges
            if length(edge_points) >= 2
                edge_indices = Vector{Int}()
                for point in edge_points
                    push!(all_points, point)
                    push!(edge_indices, point_count)
                    point_count += 1
                end
                push!(polylines, edge_indices)
                push!(polyline_bc_ids, bc_id)
            end
        end
    end

    open(vtkfile, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "BC Edges step=$step, time=$(step*Δt)")
        println(io, "ASCII")
        println(io, "DATASET POLYDATA")

        # Points
        println(io, "POINTS $(length(all_points)) float")
        for point in all_points
            if length(point) == 2
                println(io, "$(point[1]) $(point[2]) 0.0")
            else
                println(io, "$(point[1]) $(point[2]) $(point[3])")
            end
        end

        # Lines
        total_size = sum(length(poly) + 1 for poly in polylines)
        println(io, "LINES $(length(polylines)) $total_size")
        for poly in polylines
            print(io, "$(length(poly))")
            for idx in poly
                print(io, " $idx")
            end
            println(io)
        end

        # Cell data
        println(io, "CELL_DATA $(length(polylines))")
        println(io, "SCALARS bc_id int 1")
        println(io, "LOOKUP_TABLE default")
        for bc_id in polyline_bc_ids
            println(io, bc_id)
        end
    end

    pvdfile = joinpath(root, pvd_name)
    current_time = step * Δt
    
    if step == 0 || !isfile(pvdfile)
        open(pvdfile, "w") do io
            println(io, """<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">""")
            println(io, "  <Collection>")
            println(io, "    <DataSet timestep=\"" *
                        string(round(current_time, digits=6)) *
                        "\" group=\"\" part=\"0\" file=\"" *
                        basename(vtkfile) * "\"/>")
            println(io, "  </Collection>")
            println(io, "</VTKFile>")
        end
    else
        lines = readlines(pvdfile)
        insert!(lines, length(lines)-1,
                "    <DataSet timestep=\"" *
                string(round(current_time, digits=6)) *
                "\" group=\"\" part=\"0\" file=\"" *
                basename(vtkfile) * "\"/>")
        open(pvdfile, "w") do io
            for l in lines
                println(io, l)
            end
        end
    end

    return vtkfile
end

function save_bc_edges_vtk(geo::Post3D, root::String, step::Int, Δt::Float64; pvd_name="bc_surfaces.pvd")

    mkpath(root)
    new_bc_edges = geo.bc_surface
    vtkfile =  joinpath(root, "bc_surfaces_" * lpad(string(step), 4, '0') * ".vtk")

    all_points = Vector{Vector{Float64}}()
    polylines = Vector{Vector{Int}}()
    polyline_bc_ids = Vector{Int}()

    point_count = 0
    for (bc_id, edges) in new_bc_edges
        for edge_points in edges
            if length(edge_points) >= 2
                edge_indices = Vector{Int}()
                for point in edge_points
                    push!(all_points, point)
                    push!(edge_indices, point_count)
                    point_count += 1
                end
                push!(polylines, edge_indices)
                push!(polyline_bc_ids, bc_id)
            end
        end
    end

    open(vtkfile, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "BC Edges step=$step, time=$(step*Δt)")
        println(io, "ASCII")
        println(io, "DATASET POLYDATA")

        # Points
        println(io, "POINTS $(length(all_points)) float")
        for point in all_points
            if length(point) == 2
                println(io, "$(point[1]) $(point[2]) 0.0")
            else
                println(io, "$(point[1]) $(point[2]) $(point[3])")
            end
        end

        # Lines
        total_size = sum(length(poly) + 1 for poly in polylines)
        println(io, "LINES $(length(polylines)) $total_size")
        for poly in polylines
            print(io, "$(length(poly))")
            for idx in poly
                print(io, " $idx")
            end
            println(io)
        end

        # Cell data
        println(io, "CELL_DATA $(length(polylines))")
        println(io, "SCALARS bc_id int 1")
        println(io, "LOOKUP_TABLE default")
        for bc_id in polyline_bc_ids
            println(io, bc_id)
        end
    end

    pvdfile = joinpath(root, pvd_name)
    current_time = step * Δt
    
    if step == 0 || !isfile(pvdfile)
        open(pvdfile, "w") do io
            println(io, """<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">""")
            println(io, "  <Collection>")
            println(io, "    <DataSet timestep=\"" *
                        string(round(current_time, digits=6)) *
                        "\" group=\"\" part=\"0\" file=\"" *
                        basename(vtkfile) * "\"/>")
            println(io, "  </Collection>")
            println(io, "</VTKFile>")
        end
    else
        lines = readlines(pvdfile)
        insert!(lines, length(lines)-1,
                "    <DataSet timestep=\"" *
                string(round(current_time, digits=6)) *
                "\" group=\"\" part=\"0\" file=\"" *
                basename(vtkfile) * "\"/>")
        open(pvdfile, "w") do io
            for l in lines
                println(io, l)
            end
        end
    end

    return vtkfile
end

