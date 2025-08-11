const JOB_KWARGS = (:path, :freq, :fields)

struct FSI_job{S<:Peridynamics.AbstractSpatialSetup, FT<:AbstractFlowTimeSolver, ST<:Peridynamics.AbstractTimeSolver, O<:Peridynamics.AbstractJobOptions}
    flow_setup::KitBase.SolverSet
    spatial_setup::S
    f_time_solver::FT 
    s_time_solver::ST 
    options::O

    function FSI_job(flow_setup::KitBase.SolverSet, spatial_setup::S, f_time_solver::FT, s_time_solver::ST, options::O) where {S, FT, ST, O}
        Peridynamics.pre_submission_check(spatial_setup)
        return new{S, FT, ST, O}(flow_setup, spatial_setup, f_time_solver, s_time_solver, options)
    end
end

function FSI_job(flow_setup::KitBase.SolverSet, spatial_setup::S, f_time_solver::FT, s_time_solver::ST; kwargs...) where {S, FT, ST}
    o = Dict{Symbol,Any}(kwargs)
    Peridynamics.check_kwargs(o, JOB_KWARGS)
    options = fsi_get_job_options(spatial_setup, s_time_solver, o)
    return FSI_job(flow_setup, spatial_setup, f_time_solver, s_time_solver, options)
end

function Base.show(io::IO, @nospecialize(job::FSI_job))
    n_flow_points = Peridynamics.n_points(job.spatial_setup)  ##### in ks   
    n_solid_points = Peridynamics.n_points(job.spatial_setup)
    if job.spatial_setup isa AbstractMultibodySetup
        job_descr = "-point multibody Job with "
    else
        job_descr = "-point Job with "
    end
    solver = typeof(job.time_solver)
    print(io, n_flow_points, n_solid_points, job_descr, solver, " solver")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(job::FSI_job))
    if get(io, :compact, false)
        show(io, job)
    else
        println(io, "Job:")
        print(io, Peridynamics.msg_fields(job))
    end
    return nothing
end

function fsi_log_timesolver(options::Peridynamics.AbstractJobOptions, fv::Flowstep, vv::Thermstep)
    msg = "COUPLING SOLVER\n"
    msg *= Peridynamics.msg_qty("number of time steps", max(fv.n_steps, vv.n_steps))
    msg *= Peridynamics.msg_qty("stable time step size", min(fv.Δt, vv.Δt))
    msg *= Peridynamics.msg_qty("simulation time", min(fv.end_time, vv.end_time))    
    msg *= "FUILD SOLVER\n"
    msg *= Peridynamics.msg_qty("number of time steps ", fv.n_steps)
    msg *= Peridynamics.msg_qty("stable time step size", fv.Δt)
    msg *= Peridynamics.msg_qty("time step safety factor", fv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", fv.end_time)
    msg *= "VELOCITY VERLET TIME SOLVERS FOR TRUCTURE\n"
    msg *= Peridynamics.msg_qty("number of time steps", vv.n_steps)
    msg *= Peridynamics.msg_qty("stable time step size", vv.Δt)
    msg *= Peridynamics.msg_qty("time step safety factor", vv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", vv.end_time)
    Peridynamics.log_it(options, msg)
    return nothing
end

function fsi_log_timesolver(options::Peridynamics.AbstractJobOptions, fv::Flowstep, vv::Thermomechstep)
    msg = "COUPLING SOLVER\n"
    msg *= Peridynamics.msg_qty("number of time steps", max(fv.n_steps, vv.n_steps))
    msg *= Peridynamics.msg_qty("stable time step size", min(fv.Δt, vv.Δt))
    msg *= Peridynamics.msg_qty("simulation time", min(fv.end_time, vv.end_time))    
    msg *= "FUILD SOLVER\n"
    msg *= Peridynamics.msg_qty("number of time steps", fv.n_steps)
    msg *= Peridynamics.msg_qty("stable time step size", fv.Δt)
    msg *= Peridynamics.msg_qty("time step safety factor", fv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", fv.end_time)
    msg *= "THERMOMECH TIME SOLVERS FOR TRUCTURE\n"
    msg *= Peridynamics.msg_qty("number of time steps", vv.n_steps)
    msg *= Peridynamics.msg_qty("stable time step size for coupling_solvers", vv.Δt)
    msg *= Peridynamics.msg_qty("stable time step size for motion", vv.Δt_mech)
    msg *= Peridynamics.msg_qty("stable time step size for thermal diffusion", vv.Δt_therm)
    msg *= Peridynamics.msg_qty("time step safety factor", vv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", vv.end_time)
    Peridynamics.log_it(options, msg)
    return nothing
end

function fsi_log_timesolver(options::Peridynamics.AbstractJobOptions, fv::Flowstep, vv::Dualstep)
    msg = "COUPLING SOLVER\n"
    msg *= Peridynamics.msg_qty("number of time steps", max(fv.n_steps, vv.n_steps))
    msg *= Peridynamics.msg_qty("stable time step size", min(fv.Δt, vv.Δt))
    msg *= Peridynamics.msg_qty("simulation time", min(fv.end_time, vv.end_time))    
    msg *= "FUILD SOLVER\n"
    msg *= Peridynamics.msg_qty("number of time steps", fv.n_steps)
    msg *= Peridynamics.msg_qty("stable time step size", fv.Δt)
    msg *= Peridynamics.msg_qty("time step safety factor", fv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", fv.end_time)
    msg *= "THERMOMECH DUAL TIME SOLVERS FOR TRUCTURE\n"
    msg *= Peridynamics.msg_qty("number of time steps", vv.n_steps)
    msg *= Peridynamics.msg_qty("stable time step size", vv.Δt)
    msg *= Peridynamics.msg_qty("time step safety factor", vv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", vv.end_time)
    msg *= Peridynamics.msg_qty("maximum steps of ADR", vv.ADRn_steps)
    msg *= Peridynamics.msg_qty("Damping factor of ADR", vv.Λ)
    msg *= Peridynamics.msg_qty("Tolerance of ADR", vv.ADRerror)
    Peridynamics.log_it(options, msg)
    return nothing
end
mutable struct Bcstruct
    hsource::Matrix   # q/V on the solid wall
    pressure::Matrix  # p/V on the solid wall
    function Bcstruct(position::AbstractMatrix)
        total_points = size(position, 2)
        hsource = zeros(1, total_points)
        pressure = zeros(3, total_points)
        new(hsource, pressure)
    end
end

mutable struct Refactor
    tem_ref::Float64
    rho_ref::Float64
    u_ref::Float64
    l_ref::Float64
    p_ref::Float64
    t_ref::Float64
    
    function Refactor(t)
    
    end
end

function flow_log_spatial_setup(options::Peridynamics.AbstractJobOptions, fs::KitBase.SolverSet)
    msg = "FLOW SPATIAL SETUP\n"
    msg *= Peridynamics.msg_qty("number of flow points", fs.ps.nx*fs.ps.ny)
    minmax_x = Peridynamics.@sprintf("%.7g, %.7g", fs.ps.x0, fs.ps.x1)
    minmax_y = Peridynamics.@sprintf("%.7g, %.7g", fs.ps.y0, fs.ps.y1)
    msg *= Peridynamics.msg_qty("min, max values x-direction", minmax_x; indentation=4)
    msg *= Peridynamics.msg_qty("min, max values y-direction", minmax_y; indentation=4)
    msg *= "  INITIAL GAS PROPERTIES\n"
    msg *= Peridynamics.msg_qty("Knudsen Number (Kn)", fs.gas.Kn; indentation=4 )
    msg *= Peridynamics.msg_qty("Mach Number (Ma)", fs.gas.Ma; indentation=4)
    msg *= Peridynamics.msg_qty("Reynolds Number (Re)", fs.gas.Ma * KitBase.sound_speed(1.0, fs.gas.γ) / fs.gas.μᵣ; indentation=4)
    msg *= "  FLOW BOUNDARY CONDITIONS\n"
    for (i, side) in enumerate(["Left", "Right", "Top", "Bottom"])
        msg *= Peridynamics.msg_qty("  $side", fs.set.boundary[i]; indentation=4)
    end
    Peridynamics.log_it(options, msg)

end

struct JobOptions_fsi{F,V} <: Peridynamics.AbstractJobOptions
    export_allowed::Bool
    root::String
    vtk::String
    logfile::String
    freq::Int
    fields::F
    vtk_filebase::V
    flow::String
end
    
function fsi_get_job_options(spatial_setup, solver, o)
    root, freq = Peridynamics.get_root_and_freq(o)

    fields = Peridynamics.get_export_fields(spatial_setup, solver, o)
    vtk_filebase = get_vtk_filebase_fsi(spatial_setup, root)
    flow = joinpath(root, "flow")

    if isempty(root)
        options = Peridynamics.JobOptions(spatial_setup)
    else
        options = JobOptions_fsi(root, freq, fields, vtk_filebase, flow)
    end

    return options
end

function get_vtk_filebase_fsi(body::Peridynamics.AbstractBody, root::AbstractString)
    body_name = replace(string(Peridynamics.get_name(body)), " " => "_")
    filebase = isempty(body_name) ? "timestep" : body_name * "_timestep"
    vtk_filebase::String = joinpath(root, "structure", filebase)
    return vtk_filebase
end

function JobOptions_fsi(root::String, freq::Int, fields, vtk_filebase, flow::String)
    vtk = joinpath(root, "structure")
    flow =  joinpath(root, "flow")
    logfile = joinpath(root, "logfile.log")
    return JobOptions_fsi(true, root, vtk, logfile, freq, fields, vtk_filebase, flow)
end