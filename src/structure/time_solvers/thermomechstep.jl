"""
    Thermomechstep(; kwargs...) for theemom_timestep, MPI is yet not introduced.
"""
mutable struct Thermomechstep <: Peridynamics.AbstractTimeSolver
    end_time::Float64
    n_steps::Int
    Δt::Float64
    safety_factor::Float64
    Δt_therm::Float64
    Δt_mech::Float64

    function Thermomechstep(; time::Real=-1, steps::Int=-1, stepsize::Real=-1, safety_factor::Real=0.7,
              stepsize_t::Real=-1, stepsize_m::Real=-1)
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
         
        new(time, steps, stepsize, safety_factor, stepsize_t, stepsize_m)
    end
end

function Base.show(io::IO, @nospecialize(vv::Thermomechstep))
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

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(vv::Thermomechstep))
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

function Peridynamics.init_time_solver!(vv::Thermomechstep, dh::Peridynamics.AbstractDataHandler)
    vv.Δt_therm = calc_stable_timestep_th(dh, vv.safety_factor)
    vv.Δt_mech = Peridynamics.calc_stable_timestep(dh, vv.safety_factor)
    if vv.Δt < 0
        vv.Δt = minimum([vv.Δt_therm, vv.Δt_mech])
    elseif vv.Δt > minimum([vv.Δt_therm, vv.Δt_mech]) 
        vv.Δt =  minimum([vv.Δt_therm, vv.Δt_mech])
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
    Thermomechstep_check(vv)
    return nothing
end

function Thermomechstep_check(vv::Thermomechstep)
    if vv.end_time < 0
        error("`end_time` of Thermomechstep smaller than zero!\n")
    end
    if vv.n_steps < 0
        error("`n_steps` of Thermomechstep smaller than zero!\n")
    end
    if vv.Δt < 0
        error("`Δt` of Thermomechstep smaller than zero!\n")
    end
    if !(0 < vv.safety_factor < 1)
        error("`safety_factor` of Thermomechstep has invalid value!\n")
    end

    return nothing
end

function solve_thermomech_struct!(dh::Peridynamics.AbstractDataHandler, job::Job)
    options = job.options
    Δt = job.time_solver.Δt

    t_factors, m_factors = new_modify_tm(dh)
    conv, radi = find_sec_bcs_points(dh)    
    
    if mpi_isroot()
        pro = Progress(job.time_solver.n_steps; dt=1, desc="solve...", color=:normal, barlen=20,
                 enabled=Peridynamics.progress_bars())    
    end

    for n in 1:job.time_solver.n_steps
        thermomechstep_pd!(dh, options, Δt, n, t_factors, m_factors, conv, radi)
        mpi_isroot() && next!(pro)     
    end
    mpi_isroot() && finish!(pro)
    return dh
end

function thermomechstep_pd!(dh::Peridynamics.AbstractThreadsBodyDataHandler, options::Peridynamics.AbstractJobOptions,
                          Δt::Float64, n::Int, t_factors::Vector{Vector{Float64}}, m_factors::Vector{Vector{Float64}},
                          conv::Vector{Vector{Int}}, radi::Vector{Vector{Int}})
    t = n * Δt
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        Peridynamics.apply_boundary_conditions!(chunk, t)
        second_bcs!(chunk, conv[chunk_id], radi[chunk_id])
    end
    
    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_loc_to_halo!(dh, chunk_id)
        cou_calc_pflux!(dh.chunks[chunk_id], t_factors[chunk_id]) 
        new_calc_force_density!(dh.chunks[chunk_id], m_factors[chunk_id], t, Δt)
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_halo_to_loc!(dh, chunk_id)
        chunk = dh.chunks[chunk_id]        
        update_temperature!(chunk, Δt)
        Peridynamics.update_acc_and_vel!(chunk, 0.5*Δt)
        Peridynamics.update_vel_half!(chunk, 0.5*Δt)
        Peridynamics.update_disp_and_pos!(chunk, Δt)
        Peridynamics.apply_boundary_conditions!(chunk, t)
        Peridynamics.export_results(dh, options, chunk_id, n, t)
    end

    return nothing
end

function Peridynamics.req_point_data_fields_timesolver(::Type{Thermomechstep})
    fields = (:position, :displacement, :velocity, :velocity_half,
            :b_int, :b_ext, :temperature, :pflux, :hsource)
    return fields
end

function Peridynamics.req_data_fields_timesolver(::Type{Thermomechstep})
    return ()
end

function Peridynamics.log_timesolver(options::Peridynamics.AbstractJobOptions, vv::Thermomechstep)
    msg = "Thermomechstep TIME SOLVER\n"
    msg *= Peridynamics.msg_qty("number of time steps", vv.n_steps)
    msg *= Peridynamics.msg_qty("time step size", vv.Δt)
    msg *= Peridynamics.msg_qty("time step safety factor", vv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", vv.end_time)
    Peridynamics.log_it(options, msg)
    return nothing
end

