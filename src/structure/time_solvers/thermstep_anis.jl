"""
    Thermstep_anis(; kwargs...) for Anisotropic_thermal_timestep.
"""
mutable struct Thermstep_anis <: Peridynamics.AbstractTimeSolver
    end_time::Float64
    n_steps::Int
    Δt::Float64
    safety_factor::Float64

    function Thermstep_anis(; time::Real=-1, steps::Int=-1, stepsize::Real=-1, safety_factor::Real=0.7)
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
        new(time, steps, stepsize, safety_factor)
    end
end

function Base.show(io::IO, @nospecialize(vv::Thermstep_anis))
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

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(vv::Thermstep_anis))
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

function Peridynamics.init_time_solver!(vv::Thermstep_anis, dh::Peridynamics.AbstractDataHandler)
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

function Thermstep_check(vv::Thermstep_anis)
    if vv.end_time < 0
        error("`end_time` of Thermstep smaller than zero!\n")
    end
    if vv.n_steps < 0
        error("`n_steps` of Thermstep smaller than zero!\n")
    end
    if vv.Δt < 0
        error("`Δt` of Thermstep smaller than zero!\n")
    end
    return nothing
end


function solve_therm_struct_anis!(dh::Peridynamics.AbstractDataHandler, job::Job)
    options = job.options
    Δt = job.time_solver.Δt

    conv, radi = find_sec_bcs_points(dh)    
    
    if mpi_isroot()
        pro = Progress(job.time_solver.n_steps; dt=1, desc="solve...", color=:normal, barlen=20,
                 enabled=Peridynamics.progress_bars())    
    end

    for n in 1:job.time_solver.n_steps
        thermstep_pd!(dh, options, Δt, n, conv, radi)
        mpi_isroot() && next!(pro)     
    end
    mpi_isroot() && finish!(pro)
    return dh
end

function thermstep_pd!(dh::Peridynamics.AbstractThreadsBodyDataHandler, options::Peridynamics.AbstractJobOptions,
                          Δt::Float64, n::Int, conv::Vector{Vector{Int}}, radi::Vector{Vector{Int}})
    t = n * Δt
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        Peridynamics.apply_boundary_conditions!(chunk, t)
        second_bcs!(chunk, conv[chunk_id], radi[chunk_id])
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_loc_to_halo!(dh, chunk_id)
        calc_pflux!(dh.chunks[chunk_id]) 
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_halo_to_loc!(dh, chunk_id)
        chunk = dh.chunks[chunk_id]        
        update_temperature!(chunk, Δt)
        Peridynamics.export_results(dh, options, chunk_id, n, t)
    end
    return nothing
end

function Peridynamics.req_point_data_fields_timesolver(::Type{<:Thermstep_anis})
    fields = (:position, :temperature, :pflux, :hsource)
    return fields
end

function Peridynamics.req_bond_data_fields_timesolver(::Type{<:Thermstep_anis})
    return ()
end

function Peridynamics.req_data_fields_timesolver(::Type{<:Thermstep_anis})
    return ()
end

function Peridynamics.log_timesolver(options::Peridynamics.AbstractJobOptions, vv::Thermstep_anis)
    msg = "VELOCITY VERLET TIME SOLVER\n"
    msg *= Peridynamics.msg_qty("number of time steps", vv.n_steps)
    msg *= Peridynamics.msg_qty("time step size", vv.Δt)
    msg *= Peridynamics.msg_qty("time step safety factor", vv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", vv.end_time)
    Peridynamics.log_it(options, msg)
    return nothing
end

