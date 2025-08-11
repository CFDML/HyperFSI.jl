mutable struct Flowstep <: AbstractFlowTimeSolver
    end_time::Float64
    n_steps::Int
    Δt::Float64
    safety_factor::Float64

    function Flowstep(; time::Real=-1, steps::Int=-1, stepsize::Real=-1, safety_factor::Real=0.7)
        if time > 0 && steps > 0
            msg = "specify either time or number of steps, not both!"
            throw(ArgumentError(msg))
        elseif time < 0 && steps < 0
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

function init_flow_time_solver!(vv::Flowstep)
    if vv.Δt < 0
        msg = "stepsize must be given!"
        throw(ArgumentError(msg))        
    end
    if vv.end_time < 0
        vv.end_time = vv.n_steps * vv.Δt
    elseif vv.n_steps < 0
        vv.n_steps = vv.end_time ÷ vv.Δt + 1
    end
    Flowstep_check(vv)
    return nothing
end

function Flowstep_check(vv::Flowstep)
    if vv.end_time < 0
        error("`end_time` of Flowstep smaller than zero!\n")
    end
    if vv.n_steps < 0
        error("`n_steps` of Flowstep smaller than zero!\n")
    end
    if vv.Δt < 0
        error("`Δt` of Flowstep smaller than zero!\n")
    end
    return nothing
end

