### ablation function
function update_ablation_exist!(chunk::Peridynamics.AbstractBodyChunk, Δt::Float64)
    new_ablation_point_idx = Vector{Int}()
    old_ablation_point_idx = Vector{Int}()
    for i in eachindex(chunk.system.chunk_handler.loc_points)
        if chunk.storage.exist[1, i] == 1
            Δv = cal_ablation_rate(chunk, i)
            chunk.storage.ablation_deep[1, i] += Δv * Δt
    
            if chunk.storage.ablation_deep[1, i] >= chunk.system.volume[i]
                chunk.storage.exist[1, i] = 0
                id_ch = chunk.system.chunk_handler.loc_points[i] # global index
                push!(new_ablation_point_idx, id_ch)
            end
        else
            id_ch = chunk.system.chunk_handler.loc_points[i] # global index
            push!(old_ablation_point_idx, id_ch)

        end
    end
    return new_ablation_point_idx, old_ablation_point_idx
end

function cal_ablation_rate(chunk::Peridynamics.AbstractBodyChunk, i::Int) 
     (; paramsetup, storage) = chunk
    if paramsetup isa Peridynamics.AbstractPointParameters
        ablation_v = paramsetup.ablation_v
    else
        params_i = Peridynamics.get_params(paramsetup, i)
        ablation_v = params_i.ablation_v
    end

    Δv = ablation_v(storage.temperature[1, i])

    return Δv
end
#=
function update_ab_bcs_add!(chunk::Peridynamics.AbstractBodyChunk, new_ablation_point_idx::Vector{Int})
    if !isempty(new_ablation_point_idx)
        #position = chunk.system.position[:, 1:chunk.system.chunk_handler.n_loc_points]
        position = chunk.storage.position[:, 1:chunk.system.chunk_handler.n_loc_points]
        n_radius = 1.1 * chunk.paramsetup.δ
        new_add_bcs_idx = find_neighbors_indices(position, new_ablation_point_idx, n_radius)
        q_const = 1e8

        for id in new_add_bcs_idx
            chunk.storage.hsource[1, id] = q_const * chunk.storage.exist[1, id]
            #id_ch = chunk.system.chunk_handler.loc_points[id]
            #q_matrix[1, id_ch] = q_const * chunk.storage.exist[1, id]
        end
        for id in new_ablation_point_idx
            chunk.storage.hsource[1, id] = q_const * chunk.storage.exist[1, id]
            chunk.storage.velocity_half[1, id] = 0.0
            chunk.storage.velocity_half[2, id] = 0.0

            #id_ch = chunk.system.chunk_handler.loc_points[id]
            #q_matrix[1, id_ch] = q_const * chunk.storage.exist[1, id]
        end

    end
end
=#

function find_new_bcs_idx(dh::Peridynamics.AbstractThreadsBodyDataHandler, new_ablation_point_idx_all::Vector{Vector{Int}}, 
                            old_ablation_point_idx_all::Vector{Vector{Int}}, pos::Matrix{Float64})

    index_new_ablation = collect(unique(Iterators.flatten(new_ablation_point_idx_all)))
    index_old_ablation = collect(unique(Iterators.flatten(old_ablation_point_idx_all)))
    new_add_bcs_idx = Int[]
    idx_map = Dict{Int, Vector{Int}}() # map from ablation point to its neighbors
    if !isempty(index_new_ablation)
        n_radius = 1.1 * dh.chunks[1].paramsetup.δ
        new_add_bcs_idx0, idx_map = find_neighbors_indices(pos, index_new_ablation, n_radius)
        new_add_bcs_idx = setdiff(new_add_bcs_idx0, index_old_ablation) # remove the old ablation points
    end
    return new_add_bcs_idx, idx_map
end

# new boundary occurs after ablation, energy from ablation is divided by neighbors!
function update_new_bcs!(chunk::Peridynamics.AbstractBodyChunk, idx_map::Dict{Int, Vector{Int}}, new_ablation_point_idx::Vector{Int})

    for id in new_ablation_point_idx 
        ids_new_bcs = idx_map[id]
        energy_add_aver = 0.0
        if length(ids_new_bcs) > 0
            energy_add_aver = chunk.paramsetup.energy 

            for j in intersect(ids_new_bcs, chunk.system.chunk_handler.loc_points)
                j_ch = chunk.system.chunk_handler.localizer[j]
                chunk.storage.hsource[1, j_ch] = energy_add_aver
            end
        end 
    end
end

function update_new_bcs!(chunk::Peridynamics.AbstractBodyChunk, new_add_bcs_idx::Vector{Int})

    for id in intersect(new_add_bcs_idx, chunk.system.chunk_handler.loc_points) 
        i_ch = chunk.system.chunk_handler.localizer[id]
        chunk.storage.hsource[1, i_ch] = chunk.paramsetup.energy
    end
end
#=
function find_neighbors_indices(position::Matrix{Float64}, new_ablation_point_idx::Vector{Int}, radius::Float64)
    nhs = GridNeighborhoodSearch{3}(search_radius=radius, n_points=size(position, 2))
    initialize_grid!(nhs, position)

    all_neighbors_indices = Int[]

    for point_id in new_ablation_point_idx
        foreach_neighbor(position, position, nhs, point_id; search_radius=radius) do i, j, _, L
            if i != j
                push!(all_neighbors_indices, j)
            end
        end
    end
    return all_neighbors_indices
end
=#

function find_neighbors_indices(position::Matrix{Float64}, new_ablation_point_idx::Vector{Int}, radius::Float64)
    nhs = GridNeighborhoodSearch{3}(search_radius=radius, n_points=size(position, 2))
    initialize_grid!(nhs, position)

    all_neighbors_set = Set{Int}() # All neighbor indices
    neighbor_map = Dict{Int, Vector{Int}}() # neighbor lists for each ablation point

    for point_id in new_ablation_point_idx
        current_neighbors = Int[] 
        foreach_neighbor(position, position, nhs, point_id; search_radius=radius) do i, j, _, L
            if i != j
                push!(all_neighbors_set, j) 
                push!(current_neighbors, j)   
            end
        end
        neighbor_map[point_id] = current_neighbors 
    end

    all_neighbors_indices = collect(all_neighbors_set)
    return all_neighbors_indices, neighbor_map
end

function update_mech_state!(chunk::Peridynamics.AbstractBodyChunk)
    for i in eachindex(chunk.system.chunk_handler.loc_points)
        if chunk.storage.exist[1, i] < 1
            chunk.storage.velocity[:, i] .= 0.0
            chunk.storage.velocity_half[:, i] .= 0.0
            chunk.storage.velocity_half_old[:, i] .= 0.0
            chunk.storage.displacement[:, i] .= 0.0
        end
    end
end


