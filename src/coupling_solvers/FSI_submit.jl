const SUBMIT_KWARGS = (:quiet,)

# Submit for FSI simulation
function FSI_submit(job::FSI_job, bcst::Bcstruct, geo::Post2D, mode::String; output=nothing, output_vars = nothing, ref = nothing, kwargs...)
    o = Dict{Symbol,Any}(kwargs)
    Peridynamics.check_kwargs(o, SUBMIT_KWARGS)
    quiet = Peridynamics.get_submit_options(o)
    Peridynamics.set_quiet!(quiet)
    
    if mode == "T" # thermal diffusion in structure
        if Peridynamics.mpi_run()
            ret = submit_mpi_therm(job, bcst.hsource, geo; output = output, output_vars = output_vars, ref = ref)
        else 
            ret = submit_threads_therm(job, bcst.hsource, geo, nthreads(); output = output, output_vars = output_vars, ref = ref)
        end
        return ret

    elseif mode == "TM" # thermomechanics in structure
        if Peridynamics.mpi_run()
            ret = submit_mpi_thermomech(job)
        else 
            ret = submit_threads_thermomech(job, bcst, geo, nthreads(); output = output, output_vars = output_vars, ref = ref)
        end
        return ret

    elseif mode == "DTM" # Dualstep_thermomechanics in structure
        if Peridynamics.mpi_run()
            ret = submit_mpi_thermomech_dual(job)
        else 
            ret = submit_threads_thermomech_dual(job, bcst, geo, nthreads(); output = output, output_vars = output_vars, ref = ref)
        end
        return ret
    end
end

# Submit for thermal diffusion in FSI
function submit_threads_therm(job::FSI_job, hsource_bc::Matrix{Float64}, geo::Post2D, n_chunks::Int; output = nothing, output_vars = nothing, ref = nothing)

    simulation_duration = @elapsed begin
        logo_init_logs(job.options)
        flow_log_spatial_setup(job.options, job.flow_setup)
        Peridynamics.log_spatial_setup(job.options, job.spatial_setup)
        Peridynamics.log_create_data_handler_start()
        dh = Peridynamics.threads_data_handler(job.spatial_setup, job.s_time_solver, n_chunks)
        init_flow_time_solver!(job.f_time_solver)
        Peridynamics.init_time_solver!(job.s_time_solver, dh)
        Peridynamics.initialize!(dh, job.s_time_solver)
        Peridynamics.log_create_data_handler_end()
        Peridynamics.log_data_handler(job.options, dh)
        fsi_log_timesolver(job.options, job.f_time_solver, job.s_time_solver)
        solve_therm!(dh, job, geo, job.options, hsource_bc; output = output, output_vars = output_vars, ref = ref)
    end
    Peridynamics.log_simulation_duration(job.options, simulation_duration)
    return dh
end

function submit_mpi_therm(job::FSI_job, hsource_bc::Matrix{Float64}, geo::Post2D; output = nothing, output_vars = nothing, ref = nothing)

    timeit_debug_enabled() && reset_timer!(TO)
    simulation_duration = @elapsed begin
        @timeit_debug TO "initialization" begin
            logo_init_logs(job.options)
            flow_log_spatial_setup(job.options, job.flow_setup)
            log_spatial_setup(job.options, job.spatial_setup)
            log_create_data_handler_start()
            dh = mpi_data_handler(job.spatial_setup, job.s_time_solver)
            init_flow_time_solver!(job.f_time_solver)
            init_time_solver!(job.s_time_solver, dh)
            initialize!(dh, job.s_time_solver)
            log_create_data_handler_end()
            log_data_handler(job.options, dh)
            fsi_log_timesolver(job.options, job.f_time_solver, job.s_time_solver)
        end

        @timeit_debug TO "solve!" begin
            solve_therm!(dh, job, geo, job.options, hsource_bc; output = output, output_vars = output_vars, ref = ref)
        end
    end
    log_simulation_duration(job.options, simulation_duration)
    if timeit_debug_enabled()
        TimerOutputs.complement!(TO)
        log_mpi_timers(job.options)
    end
end

function solve_therm!(dh::Peridynamics.AbstractDataHandler, job::FSI_job, geo::Post2D, options::Peridynamics.AbstractJobOptions, hsource_bc::Matrix{Float64}; output = nothing, output_vars = nothing, ref = nothing)
    ks = job.flow_setup
    Δt = min(job.s_time_solver.Δt, job.f_time_solver.Δt)
    sys_n = job.s_time_solver.n_steps
    
    ctr, a1face, a2face = init_fvm(job.flow_setup; structarray = true) 

    ps = job.flow_setup.ps
    if isnothing(ref)
        ib = IBM2D(ps, geo.pos, geo.bc_edge, ReferenceVariables(;length = 1.0, temperature = 1.0, velocity = 1.0, density = 1.0, pressure = 1.0))
    else
        ib = IBM2D(ps, geo.pos, geo.bc_edge, ref)
    end

    t_factors = modify_t(dh)
    conv, radi = find_sec_bcs_points(dh)    
    
    pro = Progress(sys_n; dt=1, desc="solve...", color=:normal, barlen=20,
                 enabled=Peridynamics.progress_bars())    
    for n in 1:sys_n
        cur_t = n * Δt   
        thermstep_pd!(dh, options, Δt, n, t_factors, conv, radi)

        pd_out = merge_chunk(dh)
        temperature = pd_out[size(pd_out, 1)-1,:] .+ job.spatial_setup.point_params[1].rft #[pd_out[5, hull_indices[iter]] for iter in bi2pd]  
        update_ghost!(ctr, ps, ks.gas, ib, temperature)

        ibm_step!(job.flow_setup, ib, ctr, a1face, a2face, n, Δt/ib.rv.time, options)

        update_hsource_on_edge!(ps, ctr, geo.pos, geo.area, temperature, hsource_bc, ib)

        output_to_main!(Base.@locals, options, n, output, output_vars)

        next!(pro)     
    end
    finish!(pro)

end

# Submit for thermomechanics in FSI 
function submit_threads_thermomech(job::FSI_job, bcst::Bcstruct, geo::Post2D, n_chunks::Int; output = nothing, output_vars = nothing, ref = nothing)

    simulation_duration = @elapsed begin
        logo_init_logs(job.options)
        flow_log_spatial_setup(job.options, job.flow_setup)
        Peridynamics.log_spatial_setup(job.options, job.spatial_setup)
        Peridynamics.log_create_data_handler_start()
        dh = Peridynamics.threads_data_handler(job.spatial_setup, job.s_time_solver, n_chunks)
        init_flow_time_solver!(job.f_time_solver)
        Peridynamics.init_time_solver!(job.s_time_solver, dh)
        Peridynamics.initialize!(dh, job.s_time_solver)
        Peridynamics.log_create_data_handler_end()
        Peridynamics.log_data_handler(job.options, dh)
        fsi_log_timesolver(job.options, job.f_time_solver, job.s_time_solver)
        solve_thermomech!(dh, job, geo, bcst; output = output, output_vars = output_vars, ref = ref)
    end
    Peridynamics.log_simulation_duration(job.options, simulation_duration)
    return dh
end

function solve_thermomech!(dh::Peridynamics.AbstractDataHandler, job::FSI_job, geo::Post2D, bcst::Bcstruct; output = nothing, output_vars = nothing, ref = nothing)
    ks = job.flow_setup
    options = job.options
    Δt = min(job.s_time_solver.Δt, job.f_time_solver.Δt)
    sys_n = job.s_time_solver.n_steps
    hsource_bc = bcst.hsource
    pressure = bcst.pressure

    ctr, a1face, a2face = init_fvm(job.flow_setup; structarray = true) 

    ps = job.flow_setup.ps
    if isnothing(ref)
        ref = ReferenceVariables(;length = 1.0, temperature = 1.0, velocity = 1.0, density = 1.0, pressure = 1.0)       
    end
    ib0 = IBM2D(ps, geo.pos, geo.bc_edge, ref)

    t_factors, m_factors = new_modify_tm(dh)
    conv, radi = find_sec_bcs_points(dh)    
    
    pro = Progress(sys_n; dt=1, desc="solve...", color=:normal, barlen=20,
                 enabled=Peridynamics.progress_bars())    
    for n in 1:sys_n
        cur_t = n * Δt   
        thermomechstep_pd!(dh, options, Δt, n, t_factors, m_factors, conv, radi)

        pd_out = merge_chunk(dh)
        new_pos, new_bc_edge = update_pos(geo.pos, geo.bc_edge, pd_out)
        ib = IBM2D(ps, new_pos, new_bc_edge, ref)
        temperature = pd_out[size(pd_out, 1)-1,:] .+ job.spatial_setup.point_params[1].rft #[pd_out[5, hull_indices[iter]] for iter in bi2pd]  
        update_ghost!(ctr, ps, ks.gas, ib, temperature)

        ibm_step!(job.flow_setup, ib, ctr, a1face, a2face, n, Δt/ib.rv.time, options)

        ##### thermomechanics
        update_hsource_on_edge!(ps, ctr, geo.pos, geo.area, temperature, hsource_bc, ib)
        update_pressure_on_edge!(ps, ctr, geo.pos, geo.area, pressure, ib)

        output_to_main!(Base.@locals, options, n, output, output_vars)

        next!(pro)     
    end
    finish!(pro)

end

# Submit for thermomechanics——dualsteps in FSI 
function submit_threads_thermomech_dual(job::FSI_job, bcst::Bcstruct, geo::Post2D, n_chunks::Int; output = nothing, output_vars = nothing, ref = nothing)

    simulation_duration = @elapsed begin
        logo_init_logs(job.options)
        flow_log_spatial_setup(job.options, job.flow_setup)
        Peridynamics.log_spatial_setup(job.options, job.spatial_setup)
        Peridynamics.log_create_data_handler_start()
        dh = Peridynamics.threads_data_handler(job.spatial_setup, job.s_time_solver, n_chunks)
        init_flow_time_solver!(job.f_time_solver)
        Peridynamics.init_time_solver!(job.s_time_solver, dh)
        Peridynamics.initialize!(dh, job.s_time_solver)
        Peridynamics.log_create_data_handler_end()
        Peridynamics.log_data_handler(job.options, dh)
        fsi_log_timesolver(job.options, job.f_time_solver, job.s_time_solver)
        solve_thermomech_dual!(dh, job, geo, bcst; output = output, output_vars = output_vars, ref = ref)
    end
    Peridynamics.log_simulation_duration(job.options, simulation_duration)
    return dh
end

function solve_thermomech_dual!(dh::Peridynamics.AbstractDataHandler, job::FSI_job, geo::Post2D, bcst::Bcstruct; output = nothing, output_vars = nothing, ref = nothing)
    ks = job.flow_setup
    options = job.options
    Δt = min(job.s_time_solver.Δt, job.f_time_solver.Δt)
    sys_n = job.s_time_solver.n_steps
    hsource_bc = bcst.hsource
    pressure = bcst.pressure
    
    ctr, a1face, a2face = init_fvm(job.flow_setup; structarray = true) 

    ps = job.flow_setup.ps
    if isnothing(ref)
        ref = ReferenceVariables(;length = 1.0, temperature = 1.0, velocity = 1.0, density = 1.0, pressure = 1.0)       
    end
    ib0 = IBM2D(ps, geo.pos, geo.bc_edge, ref)

    t_factors, m_factors = new_modify_tm(dh)
    conv, radi = find_sec_bcs_points(dh)    
    
    pro = Progress(sys_n; dt=1, desc="solve...", color=:normal, barlen=20,
                 enabled=Peridynamics.progress_bars())    
    for n in 1:sys_n
        cur_t = n * Δt   
        thermomech_dualstep_pd!(dh, options, Δt, n, t_factors, m_factors, conv, radi, job.s_time_solver)

        pd_out = merge_chunk(dh)
        new_pos, new_bc_edge = update_pos(geo.pos, geo.bc_edge, pd_out)
        ib = IBM2D(ps, new_pos, new_bc_edge, ref)
        temperature = pd_out[size(pd_out, 1)-1,:] .+ job.spatial_setup.point_params[1].rft #[pd_out[5, hull_indices[iter]] for iter in bi2pd]  
        update_ghost!(ctr, ps, ks.gas, ib, temperature)

        ibm_step!(job.flow_setup, ib, ctr, a1face, a2face, n, Δt/ib.rv.time, options)
        ##### thermomechanics
        update_hsource_on_edge!(ps, ctr, geo.pos, geo.area, temperature, hsource_bc, ib)
        update_pressure_on_edge!(ps, ctr, geo.pos, geo.area, pressure, ib)

        output_to_main!(Base.@locals, options, n, output, output_vars)

        next!(pro)     
    end
    finish!(pro)

end

# Submit for PD simulation
function FSI_submit(job::Job, mode::String; kwargs...)
    o = Dict{Symbol,Any}(kwargs)
    Peridynamics.check_kwargs(o, SUBMIT_KWARGS)
    quiet = Peridynamics.get_submit_options(o)
    Peridynamics.set_quiet!(quiet)
    
    if mode == "T" # thermal diffusion in structure
        if Peridynamics.mpi_run()
            ret = T_submit_mpi(job)
        else 
            ret = T_submit_threads(job, nthreads())
        end
        return ret

    elseif mode == "TM" # thermomechanics in structure
        if Peridynamics.mpi_run()
            ret = TM_submit_mpi(job)
        else 
            ret = TM_submit_threads(job, nthreads())
        end
        return ret

    elseif mode == "DTM" # thermomechanics_dual in structure
        if Peridynamics.mpi_run()
            ret = DTM_submit_mpi(job)
        else 
            ret = DTM_submit_threads(job, nthreads())
        end
        return ret
    end
end

function T_submit_threads(job::Job, n_chunks::Int)

    simulation_duration = @elapsed begin
        logo_init_logs(job.options)
        Peridynamics.log_spatial_setup(job.options, job.spatial_setup)
        Peridynamics.log_create_data_handler_start()
        dh = Peridynamics.threads_data_handler(job.spatial_setup, job.time_solver, n_chunks)
        Peridynamics.init_time_solver!(job.time_solver, dh)
        Peridynamics.initialize!(dh, job.time_solver)
        Peridynamics.log_create_data_handler_end()
        Peridynamics.log_data_handler(job.options, dh)
        Peridynamics.log_timesolver(job.options, job.time_solver)
        solve_therm_struct!(dh, job)
    end
    Peridynamics.log_simulation_duration(job.options, simulation_duration)
    return dh
end

function TM_submit_threads(job::Job, n_chunks::Int)

    simulation_duration = @elapsed begin
        logo_init_logs(job.options)
        Peridynamics.log_spatial_setup(job.options, job.spatial_setup)
        Peridynamics.log_create_data_handler_start()
        dh = Peridynamics.threads_data_handler(job.spatial_setup, job.time_solver, n_chunks)
        Peridynamics.init_time_solver!(job.time_solver, dh)
        Peridynamics.initialize!(dh, job.time_solver)
        Peridynamics.log_create_data_handler_end()
        Peridynamics.log_data_handler(job.options, dh)
        Peridynamics.log_timesolver(job.options, job.time_solver)
        solve_thermomech_struct!(dh, job)
    end
    Peridynamics.log_simulation_duration(job.options, simulation_duration)
    return dh
end

function DTM_submit_threads(job::Job, n_chunks::Int)

    simulation_duration = @elapsed begin
        logo_init_logs(job.options)
        Peridynamics.log_spatial_setup(job.options, job.spatial_setup)
        Peridynamics.log_create_data_handler_start()
        dh = Peridynamics.threads_data_handler(job.spatial_setup, job.time_solver, n_chunks)
        Peridynamics.init_time_solver!(job.time_solver, dh)
        Peridynamics.initialize!(dh, job.time_solver)
        Peridynamics.log_create_data_handler_end()
        Peridynamics.log_data_handler(job.options, dh)
        Peridynamics.log_timesolver(job.options, job.time_solver)
        solve_dual_thermomech_struct!(dh, job)
    end
    Peridynamics.log_simulation_duration(job.options, simulation_duration)
    return dh
end

function merge_chunk(C::Peridynamics.AbstractThreadsBodyDataHandler)
    n = C.chunks[end].system.chunk_handler.loc_points[end]
    pd_out = fill(0.0, 8, n)
    for i in 1:C.n_chunks
        for j in C.chunks[i].system.chunk_handler.loc_points
            m = C.chunks[i].system.chunk_handler.localizer[j]
            pd_out[1:3, j] += C.chunks[i].storage.position[1:3, m]
            pd_out[4:6, j] += C.chunks[i].storage.velocity[1:3, m]
            pd_out[7, j] += C.chunks[i].storage.temperature[1, m]
            pd_out[8, j] += C.chunks[i].storage.damage[m]
        end
    end
    return pd_out
end
