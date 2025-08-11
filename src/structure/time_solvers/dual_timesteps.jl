"""
    Dualstep(; kwargs...) for dual_timestep, MPI is yet not introduced.
"""
mutable struct Dualstep <: Peridynamics.AbstractTimeSolver
    end_time::Float64
    n_steps::Int
    Δt::Float64
    safety_factor::Float64
    ADRn_steps::Int
    Λ::Float64
    ADRerror::Float64

    function Dualstep(; time::Real=-1, steps::Int=-1, stepsize::Real=-1, safety_factor::Real=0.7,
                        ADRsteps::Int=-1, d_factor::Real=1.0, ADRerror::Real=1e-5)
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
        new(time, steps, stepsize, safety_factor, ADRsteps, d_factor, ADRerror)
    end
end

function Base.show(io::IO, @nospecialize(vv::Dualstep))
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

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(vv::Dualstep))
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

function Peridynamics.init_time_solver!(vv::Dualstep, dh::Peridynamics.AbstractDataHandler)
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

function Dualstep_check(vv::Dualstep)
    if vv.end_time < 0
        error("`end_time` of Dualstep smaller than zero!\n")
    end
    if vv.n_steps < 0
        error("`n_steps` of Dualstep smaller than zero!\n")
    end
    if vv.Δt < 0
        error("`Δt` of Dualstep smaller than zero!\n")
    end
    if !(0 < vv.safety_factor < 1)
        error("`safety_factor` of Dualstep has invalid value!\n")
    end
    if vv.ADRn_steps < 0
        error("`n_steps` of Dualstep smaller than zero!\n")
    end
    if vv.Λ < 0
        error("`Λ` of Dualstep smaller than zero!\n")
    end
    if !(0 < vv.ADRerror < 1)
        error("`ADRerror` of Dualstep has invalid value!\n")
    end
    return nothing
end

function solve_dual_thermomech_struct!(dh::Peridynamics.AbstractDataHandler, job::Job)
    
    options = job.options
    vv = job.time_solver
    Δt = vv.Δt
    adr = Peridynamics.DynamicRelaxation(steps = vv.ADRn_steps, damping_factor=vv.Λ)
    Peridynamics.init_density_matrix!(dh, adr)
    t_factors, m_factors = new_modify_tm(dh)
    conv, radi = find_sec_bcs_points(dh)    
    
    if mpi_isroot()
        p = Progress(vv.n_steps; dt=1, desc="TIME INTEGRATION LOOP", color=:normal,
                     barlen=40, enabled=Peridynamics.progress_bars())
    end

    for n in 1:vv.n_steps
        thermomech_dualstep_pd!(dh, options, Δt, n, t_factors, m_factors, conv, radi, vv)
        Peridynamics.mpi_isroot() && next!(p)
    end
    Peridynamics.mpi_isroot() && Peridynamics.finish!(p)
    return dh

end

function thermomech_dualstep_pd!(dh::Peridynamics.AbstractThreadsBodyDataHandler, options::Peridynamics.AbstractJobOptions, Δt::Float64, n::Int,
                        t_factors::Vector{Vector{Float64}}, m_factors::Vector{Vector{Float64}},
                        conv::Vector{Vector{Int}}, radi::Vector{Vector{Int}}, vv::Dualstep)
    t = n * Δt
    aΔt = 1.0

    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        Peridynamics.apply_boundary_conditions!(chunk, t)
        second_bcs!(chunk, conv[chunk_id], radi[chunk_id])
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_loc_to_halo!(dh, chunk_id)
        cou_calc_pflux!(dh.chunks[chunk_id], t_factors[chunk_id]) 
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_halo_to_loc!(dh, chunk_id)
        chunk = dh.chunks[chunk_id]        
        update_temperature!(chunk, Δt)
    end

    for adrn in 1:vv.ADRn_steps     
        @threads :static for chunk_id in eachindex(dh.chunks)
            chunk = dh.chunks[chunk_id]
            Peridynamics.apply_boundary_conditions!(chunk, t)
            Peridynamics.update_disp_and_pos!(chunk, aΔt)
        end
        
        error = cal_error(dh.chunks)
        new_calc_force_density!(dh, m_factors, t)

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

function cal_error(chunks) #::Vector{Peridynamics.BodyChunk{}}
    errorv = 0.0
    erroru = 0.0
    #
    for i in eachindex(chunks)
        for j in eachindex(chunks[i].system.chunk_handler.loc_points)
            errorv += abs(chunks[i].storage.velocity_half[1,j]) 
            + abs(chunks[i].storage.velocity_half[2,j])
            + abs(chunks[i].storage.velocity_half[3,j])

            erroru += abs(chunks[i].storage.displacement[1,j]) 
            + abs(chunks[i].storage.displacement[2,j])
            + abs(chunks[i].storage.displacement[3,j])
        end
    end
    if erroru ≈ 0.0
        error = 1.0
    else
        error = errorv/erroru
    end
    #
    #=
    for i in eachindex(chunks)
        for j in eachindex(chunks[i].system.chunk_handler.loc_points)
            errorv += (chunks[i].storage.velocity_half[1,j])^2 
            + (chunks[i].storage.velocity_half[2,j])^2
            + (chunks[i].storage.velocity_half[3,j])^2
        end
    end

    error = sqrt(errorv / chunks[end].system.chunk_handler.loc_points[end])
    =#
    return error
end

function Peridynamics.req_point_data_fields_timesolver(::Type{Dualstep})
    fields = (:position, :displacement, :velocity, :velocity_half, :velocity_half_old,
            :b_int, :b_int_old, :b_ext, :density_matrix, :temperature, :pflux, :hsource)
    return fields
end

function Peridynamics.req_data_fields_timesolver(::Type{Dualstep})
    return ()
end

function Peridynamics.log_timesolver(options::Peridynamics.AbstractJobOptions, vv::Dualstep)
    msg = "Dualstep TIME SOLVER\n"
    msg *= Peridynamics.msg_qty("number of time steps", vv.n_steps)
    msg *= Peridynamics.msg_qty("time step size", vv.Δt)
    msg *= Peridynamics.msg_qty("time step safety factor", vv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", vv.end_time)
    msg *= Peridynamics.msg_qty("maximum steps of ADR", vv.ADRn_steps)
    msg *= Peridynamics.msg_qty("Damping factor of ADR", vv.Λ)
    msg *= Peridynamics.msg_qty("Tolerance of ADR", vv.ADRerror)
    Peridynamics.log_it(options, msg)
    return nothing
end

