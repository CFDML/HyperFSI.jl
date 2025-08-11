function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:position})
    return copy(system.position)
end

function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:displacement})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:velocity})
    return  zeros(3,size(system.position, 2))
end

function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:velocity_half})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:velocity_half_old})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function init_field_solver(::Peridynamics.AbstractTimeSolver, ::Peridynamics.AbstractSystem, ::Val{:velocity_half_old})
    return Array{Float64,2}(undef, 0, 0)
end

function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:acceleration})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:b_int})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:b_int_old})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function init_field_solver(::Peridynamics.AbstractTimeSolver, ::Peridynamics.AbstractSystem, ::Val{:b_int_old})
    return Array{Float64,2}(undef, 0, 0)
end

function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:b_ext})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:density_matrix})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function init_field_solver(::Peridynamics.AbstractTimeSolver, ::Peridynamics.AbstractSystem, ::Val{:density_matrix})
    return Array{Float64,2}(undef, 0, 0)
end

function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:temperature})
    return zeros(1, Peridynamics.size(system.position, 2))
end

function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:pflux})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep, system::Peridynamics.AbstractSystem, ::Val{:hsource})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end



