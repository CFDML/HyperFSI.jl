"""
    Thermstep_ablation, MPI is yet not introduced.
"""
mutable struct Thermstep_ablation <: Peridynamics.AbstractTimeSolver
    end_time::Float64
    n_steps::Int
    Δt::Float64
    safety_factor::Float64
    ablation_update_freq::Int

    function Thermstep_ablation(; time::Real=-1, steps::Int=-1, stepsize::Real=-1, safety_factor::Real=0.7, ablation_freq::Int=1)
        if time < 0 && steps < 0
            msg = "specify either time or number of steps!"
            throw(ArgumentError(msg))
        end
        if !(0 < safety_factor < 1)
            msg = "wrong safety factor specified! condition: 0 < safety_factor < 1"
            throw(ArgumentError(msg))
        end
        if stepsize > 0
            @warn "stepsize specified! Please be sure that the CFD-condition holds!"
        end  
        
        if ablation_freq < 0
            msg = "ablation_update_freq should be a positive integer!"
            throw(ArgumentError(msg))
        end
        new(time, steps, stepsize, safety_factor, ablation_freq)
    end
end

function Base.show(io::IO, @nospecialize(vv::Thermstep_ablation))
    print(io, typeof(vv))
    fields = Vector{Symbol}()
    for field in fieldnames(typeof(vv))
        value = Peridynamics.getfield(vv, field)
        if value > 0
            push!(fields, field)
        end
    end
    print(io, Peridynamics.msg_fields_in_brackets(vv, Tuple(fields)))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(vv::Thermstep_ablation))
    if get(io, :compact, false)
        show(io, vv)
    else
        println(io, typeof(vv), ":")
        fields = Vector{Symbol}()
        for field in fieldnames(typeof(vv))
            value = Peridynamics.getfield(vv, field)
            if value > 0
                push!(fields, field)
            end
        end
        print(io, Peridynamics.msg_fields(vv, Tuple(fields)))
    end
    return nothing
end

function Peridynamics.init_time_solver!(vv::Thermstep_ablation, dh::Peridynamics.AbstractDataHandler)
    if vv.Δt < 0
        vv.Δt = calc_stable_timestep_th(dh, vv.safety_factor)
    else
        vv.Δt = min(calc_stable_timestep_th(dh, vv.safety_factor), vv.Δt)
    end
    if vv.end_time < 0
        vv.end_time = vv.n_steps * vv.Δt
    elseif vv.n_steps < 0
        vv.n_steps = vv.end_time ÷ vv.Δt + 1
    elseif vv.n_steps >= 0
        if vv.Δt > vv.end_time / vv.n_steps
            vv.Δt =  vv.end_time / vv.n_steps
        else
            vv.n_steps =  vv.end_time ÷ vv.Δt + 1
        end
    end
    Thermstep_check(vv)
    return nothing
end

function Thermstep_check(vv::Thermstep_ablation)
    if vv.end_time < 0
        error("`end_time` of Thermstep_ablation smaller than zero!\n")
    end
    if vv.n_steps < 0
        error("`n_steps` of Thermstep_ablation smaller than zero!\n")
    end
    if vv.Δt < 0
        error("`Δt` of Thermstep_ablation smaller than zero!\n")
    end

    if vv.ablation_update_freq < 0
        error("`ablation_update_freq` of Thermstep_ablation should be a positive integer!\n")
    end

    if vv.ablation_update_freq > vv.n_steps
        @error("`ablation_update_freq` of Thermstep_ablation should be no more than total simulation steps.") 
    end
    return nothing
end

function solve_therm_struct_ablation!(dh::Peridynamics.AbstractDataHandler, job::Job)
    options = job.options
    Δt = job.time_solver.Δt
    ablation_update_freq = job.time_solver.ablation_update_freq

    conv, radi = find_sec_bcs_points(dh)
    position = reduce(hcat, [chunk.storage.position[:, 1:chunk.system.chunk_handler.n_loc_points] for chunk in dh.chunks])        
    
    if mpi_isroot()
        pro = Progress(job.time_solver.n_steps; dt=1, desc="solve...", color=:normal, barlen=20,
                 enabled=Peridynamics.progress_bars())    
    end

    for n in 1:job.time_solver.n_steps
        thermstep_ablation_pd!(dh, options, Δt, n, conv, radi, ablation_update_freq, position)
        mpi_isroot() && next!(pro)     
    end
    mpi_isroot() && finish!(pro)
    return dh
end

function thermstep_ablation_pd!(dh::Peridynamics.AbstractThreadsBodyDataHandler, options::Peridynamics.AbstractJobOptions,
                          Δt::Float64, n::Int, conv::Vector{Vector{Int}}, radi::Vector{Vector{Int}}, n_ablation::Int, pos::Matrix{Float64})
    t = n * Δt
    new_ab_idx_all = [Int[] for _ in 1:dh.n_chunks]
    old_ab_idx_all = [Int[] for _ in 1:dh.n_chunks]
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        Peridynamics.apply_boundary_conditions!(chunk, t)
        second_bcs!(chunk, conv[chunk_id], radi[chunk_id])
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_loc_to_halo!(dh, chunk_id)
        calc_pflux_ablation!(dh.chunks[chunk_id]) 
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_halo_to_loc!(dh, chunk_id)
        chunk = dh.chunks[chunk_id]        
        update_temperature!(chunk, Δt)
        Peridynamics.apply_boundary_conditions!(chunk, t)
        if mod(n, n_ablation) == 0
            new_ab_idx_all[chunk_id], old_ab_idx_all[chunk_id] = update_ablation_exist!(chunk, Δt)
        end
    end

    new_bcs_idx, idx_map = find_new_bcs_idx(dh, new_ab_idx_all, old_ab_idx_all, pos)

    @threads :static for chunk_id in eachindex(dh.chunks)
        #update_new_bcs!(dh.chunks[chunk_id], idx_map, new_ab_idx_all[chunk_id])
        update_new_bcs!(dh.chunks[chunk_id], new_bcs_idx)
        Peridynamics.export_results(dh, options, chunk_id, n, t)
    end

    return nothing
end

function Peridynamics.req_point_data_fields_timesolver(::Type{Thermstep_ablation})
    fields = (:position, :temperature, :pflux, :hsource)
    return fields
end

function Peridynamics.req_data_fields_timesolver(::Type{Thermstep_ablation})
    return ()
end

function Peridynamics.log_timesolver(options::Peridynamics.AbstractJobOptions, vv::Thermstep_ablation)
    msg = "VELOCITY VERLET TIME SOLVER FOR ABLATION\n"
    msg *= Peridynamics.msg_qty("number of time steps", vv.n_steps)
    msg *= Peridynamics.msg_qty("time step size", vv.Δt)
    msg *= Peridynamics.msg_qty("time step safety factor", vv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", vv.end_time)
    msg *= Peridynamics.msg_qty("ablation update freq", vv.ablation_update_freq)
    Peridynamics.log_it(options, msg)
    return nothing
end

