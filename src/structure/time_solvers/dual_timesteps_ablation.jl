"""
    Dualstep(; kwargs...) for dual_timestep, MPI is yet not introduced.
"""
mutable struct Dualstep_ablation <: Peridynamics.AbstractTimeSolver
    end_time::Float64
    n_steps::Int
    Δt::Float64
    safety_factor::Float64
    ADRn_steps::Int
    Λ::Float64
    ADRerror::Float64
    ablation_update_freq::Int

    function Dualstep_ablation(; time::Real=-1, steps::Int=-1, stepsize::Real=-1, safety_factor::Real=0.7,
                        ADRsteps::Int=-1, d_factor::Real=1.0, ADRerror::Real=1e-5, ablation_freq::Int=1)
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
        
        if ADRsteps ≤ 0
            msg = "`ADR_steps` should be larger than zero!\n"
            throw(ArgumentError(msg))
        end
        if d_factor ≤ 0
            msg = "`damping_factor` should be larger than zero!\n"
            throw(ArgumentError(msg))
        end  
        if ablation_freq < 0
            msg = "ablation_update_freq should be a positive integer!"
            throw(ArgumentError(msg))
        end              
        new(time, steps, stepsize, safety_factor, ADRsteps, d_factor, ADRerror, ablation_freq)
    end
end

function Base.show(io::IO, @nospecialize(vv::Dualstep_ablation))
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

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(vv::Dualstep_ablation))
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

function Peridynamics.init_time_solver!(vv::Dualstep_ablation, dh::Peridynamics.AbstractDataHandler)
    if vv.Δt < 0
        vv.Δt = calc_stable_timestep_th(dh, vv.safety_factor)
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
    Dualstep_check(vv)
    return nothing
end

function Dualstep_check(vv::Dualstep_ablation)
    if vv.end_time < 0
        error("`end_time` of Dualstep_ablation smaller than zero!\n")
    end
    if vv.n_steps < 0
        error("`n_steps` of Dualstep_ablation smaller than zero!\n")
    end
    if vv.Δt < 0
        error("`Δt` of Dualstep_ablation smaller than zero!\n")
    end
    if !(0 < vv.safety_factor < 1)
        error("`safety_factor` of Dualstep_ablation has invalid value!\n")
    end
    if vv.ADRn_steps < 0
        error("`n_steps` of Dualstep_ablation smaller than zero!\n")
    end
    if vv.Λ < 0
        error("`Λ` of Dualstep_ablation smaller than zero!\n")
    end
    if !(0 < vv.ADRerror < 1)
        error("`ADRerror` of Dualstep_ablation has invalid value!\n")
    end
    if vv.ablation_update_freq > vv.n_steps
        @error("`ablation_update_freq` of Thermstep_ablation should be no more than total simulation steps.") 
    end    
    return nothing
end

function solve_dual_thermomech_struct_ablation!(dh::Peridynamics.AbstractDataHandler, job::Job, geo::AbstractPDGeometry)
    
    options = job.options
    vv = job.time_solver
    Δt = vv.Δt
    ablation_update_freq = vv.ablation_update_freq

    adr = Peridynamics.DynamicRelaxation(steps = vv.ADRn_steps, damping_factor=vv.Λ)
    Peridynamics.init_density_matrix!(dh, adr)
    conv, radi = find_sec_bcs_points(dh)    
    
    if mpi_isroot()
        p = Progress(vv.n_steps; dt=1, desc="TIME INTEGRATION LOOP", color=:normal,
                     barlen=40, enabled=Peridynamics.progress_bars())
    end

    for n in 1:vv.n_steps
        thermomech_Dualstep_ablation_pd!(dh, options, Δt, n, conv, radi, vv, ablation_update_freq, geo)
        Peridynamics.mpi_isroot() && next!(p)
    end
    Peridynamics.mpi_isroot() && Peridynamics.finish!(p)
    return dh

end

function thermomech_Dualstep_ablation_pd!(dh::Peridynamics.AbstractThreadsBodyDataHandler, options::Peridynamics.AbstractJobOptions, Δt::Float64, n::Int,
                        conv::Vector{Vector{Int}}, radi::Vector{Vector{Int}}, vv::Dualstep_ablation, n_ablation::Int, geo::AbstractPDGeometry)
    t = n * Δt
    aΔt = 1.0 # time step size of quasi-static mechanical step 
    new_ab_idx_all = [Int[] for _ in 1:dh.n_chunks]
    old_ab_idx_all = [Int[] for _ in 1:dh.n_chunks]    

    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        Peridynamics.apply_boundary_conditions!(chunk, t)
        second_bcs!(chunk, conv[chunk_id], radi[chunk_id])
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_loc_to_halo!(dh, chunk_id)
        cou_calc_pflux!(dh.chunks[chunk_id]) 
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_halo_to_loc!(dh, chunk_id)
        chunk = dh.chunks[chunk_id]        
        update_temperature!(chunk, Δt)
        new_ab_idx_all[chunk_id], old_ab_idx_all[chunk_id] = update_ablation_exist!(chunk, Δt)        
    end

    index_new_ablation = collect(unique(Iterators.flatten(new_ab_idx_all)))
    if !isempty(index_new_ablation) && (mod(n, n_ablation) == 0)   
        index_old_ablation = collect(unique(Iterators.flatten(old_ab_idx_all)))
        new_bcs_idx, idx_map = find_new_bcs_idx(dh, index_new_ablation, index_old_ablation, geo.pos)

        @threads :static for chunk_id in eachindex(dh.chunks)
            update_new_bcs_f!(dh.chunks[chunk_id], idx_map, new_ab_idx_all[chunk_id]) #new force BCs
        end
        update_new_edges!(geo, index_new_ablation)
    end  
    if mod(n, n_ablation) == 0
        save_bc_edges_vtk(geo, options.root, n, Δt)
    end

    for adrn in 1:vv.ADRn_steps     
        @threads :static for chunk_id in eachindex(dh.chunks)
            chunk = dh.chunks[chunk_id]
            Peridynamics.apply_boundary_conditions!(chunk, t)
            Peridynamics.update_disp_and_pos!(chunk, aΔt)
        end
        
        error = cal_error(dh.chunks)
        new_calc_force_density_ablation!(dh, t, aΔt)

        @threads :static for chunk_id in eachindex(dh.chunks)
            chunk = dh.chunks[chunk_id]
            cn = Peridynamics.calc_damping(chunk, aΔt)
            if adrn == 1 #&& n == 1
                Peridynamics.relaxation_first_step!(chunk, aΔt)
            else
                Peridynamics.relaxation_step!(chunk, aΔt, cn)
            end
        end   
        
        if error <= vv.ADRerror
            println("ADR convergence :)")
            break
        end
        
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.export_results(dh, options, chunk_id, n, t)
    end

    return nothing
end

function Peridynamics.req_point_data_fields_timesolver(::Type{Dualstep_ablation})
    fields = (:position, :displacement, :velocity, :velocity_half, :velocity_half_old,
            :b_int, :b_int_old, :b_ext, :density_matrix, :temperature, :pflux, :hsource, :exist)
    return fields
end

function Peridynamics.req_data_fields_timesolver(::Type{Dualstep_ablation})
    return ()
end

function Peridynamics.log_timesolver(options::Peridynamics.AbstractJobOptions, vv::Dualstep_ablation)
    msg = "Dualstep TIME SOLVER for Ablation\n"
    msg *= Peridynamics.msg_qty("number of time steps", vv.n_steps)
    msg *= Peridynamics.msg_qty("time step size", vv.Δt)
    msg *= Peridynamics.msg_qty("time step safety factor", vv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", vv.end_time)
    msg *= Peridynamics.msg_qty("maximum steps of ADR", vv.ADRn_steps)
    msg *= Peridynamics.msg_qty("Damping factor of ADR", vv.Λ)
    msg *= Peridynamics.msg_qty("Tolerance of ADR", vv.ADRerror)
    msg *= Peridynamics.msg_qty("ablation update freq", vv.ablation_update_freq)    
    Peridynamics.log_it(options, msg)
    return nothing
end

