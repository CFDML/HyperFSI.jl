"""
Bond-Based dualstep_thermomechanics Module with Temperature-Dependent Properties and Ablation

This module extends the bond-based dualstep_thermomechanics model 
by incorporating Ablation Behavior.

··This module provides functionality for ablation, primarily for material point removal.
··Uses an 'exist' state to indicate ablation status, which relates to ablation amount.
··Specific chemical processes are configured in the ablation file, with flexible complexity.
··Core functionality identifies new boundary points, 
  and more critically, establishes new watertight boundaries for subsequent fluid-structure interaction.
··Ablated points can be modeled to describe subsequent effects.
··The variation of properties with temperature is under consideration.

··Basically, Materials model are as same as BBTMAMaterial, only times itergral is different.
"""
function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:position})
    return copy(system.position)
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:displacement})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:velocity})
    return  zeros(3,size(system.position, 2))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:velocity_half})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:velocity_half_old})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:acceleration})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:b_int})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:b_int_old})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:b_ext})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:density_matrix})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:temperature})
    return zeros(1, Peridynamics.size(system.position, 2))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:pflux})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:hsource})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:ablation_deep})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Dualstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:exist})
    return ones(Int, 1, size(system.position, 2))
end





