"""
Bond-Based Thermal Diffusion Module with Temperature-Dependent Properties and Ablation

This module extends the temperature-Dependent bond-based thermal diffusion model 
by incorporating Ablation Behavior.

··This module provides functionality for ablation, primarily for material point removal.
··Uses an 'exist' state to indicate ablation status, which relates to ablation amount.
··Specific chemical processes are configured in the ablation file, with flexible complexity.
··Core functionality identifies new boundary points, 
  and more critically, establishes new watertight boundaries for subsequent fluid-structure interaction.
··Ablated points can be modeled to describe subsequent effects.

··BBTTA means Bond-Based Thermal Diffusion with Temperature-Dependent Properties and Ablation.
"""
@inline thermal_temp_ablation_kwargs() = (:E, :nu, :G, :K, :lambda, :mu, :rho, :horizon, 
:Gc, :epsilon_c, :kc, :aph, :cv, :rft, :h, :hσ, :hϵ, :tem∞, :thick, :kp_T, :cv_T, :aph_T, :ablation_v, :energy)

struct BBTTAMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
    dmgmodel::DM
    function BBTTAMaterial{C}(dmgmodel::DM) where {C,DM} #BB_Temperature_dependent_Theermal
        new{C,DM}(dmgmodel)
    end
end

function BBTTAMaterial{C}(; dmgmodel::Peridynamics.AbstractDamageModel=CriticalStretch()) where {C}
    return BBTTAMaterial{C}(dmgmodel)
end

BBTTAMaterial(; kwargs...) = BBTTAMaterial{NoCorrection}(; kwargs...)

struct BBTTAPointParameters <: Peridynamics.AbstractPointParameters
    δ::Float64
    rho::Float64
    E::Float64
    nu::Float64
    G::Float64
    K::Float64
    λ::Float64
    μ::Float64
    Gc::Float64
    εc::Float64
    bc::Float64
    kc::Float64 # thermal conductivity at room temperature
    kp::Float64 # microconductivity at room temperature
    aph::Float64 # thermal expansion at room temperature
    cv::Float64 # specific heat capacity at room temperature
    rft::Float64 # Reference temperature at room temperature
    h::Float64 # convective heat transfer coefficient,
    hσ::Float64 # Stefan-Boltzman constant,
    hϵ::Float64 # emissivity
    tem∞::Float64 # temperature of the surrounding medium
    kp_T::Function # function describing kp variation with temperature f(kp,T)
    cv_T::Function # function describing kp variation with temperature f(kp,T)
    aph_T::Function # function describing kp variation with temperature f(kp,T)
    ablation_v::Function # function describing ablation rate f(vol_reference, T)
    energy::Float64 # total energy release from ablation
end

function const_ablation_v(T)
    return 0.0
end

function BBTTAPointParameters(mat::BBTTAMaterial, p::Dict{Symbol,Any})
    par =  Peridynamics.get_given_elastic_params(p)
    (; E, nu, G, K, λ, μ) = par

    if haskey(p, :thick)
        p[:nu] = 1/3
    else 
        p[:nu] = 1/4    
    end
   
    (; δ, rho, E, nu, G, K, λ, μ) =  Peridynamics.get_required_point_parameters(mat, p)
    Gc, εc =  Peridynamics.get_frac_params(mat.dmgmodel, p, δ, K)
    kc, aph, cv, rft, h, hσ, hϵ, tem∞, energy = get_thermal_ablation_params(p, δ) 

    if haskey(p, :thick) #2D 
        thick = float(p[:thick])
        bc = 9 * E / (π * thick * δ^3) # bond constant
        kp = 6 * kc / (π * thick * δ^3) # microcndicitvity constant
    else #3D
        bc = 12 * E / (π * δ^4) 
        kp = 6 * kc / (π * δ^4) 
    end   

    kp_T = get(p, :kp_T, const_kp)
    cv_T = get(p, :cv_T, const_cv)
    aph_T = get(p, :aph_T, const_aph)
    ablation_v = get(p, :ablation_v, const_ablation_v)

    return BBTTAPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc, kc, kp, aph, cv, rft, h, hσ, hϵ, tem∞, kp_T, cv_T, aph_T, ablation_v, energy)
end

function get_thermal_ablation_params(p::Dict{Symbol,Any}, δ)
    haskey(p, :kc) || throw(UndefKeywordError(:kc))
    haskey(p, :aph) || throw(UndefKeywordError(:aph))
    haskey(p, :cv) || throw(UndefKeywordError(:cv))
    haskey(p, :rft) || throw(UndefKeywordError(:rft))
    haskey(p, :h) || throw(UndefKeywordError(:h))  
    haskey(p, :hσ) || throw(UndefKeywordError(:hσ))  
    haskey(p, :hϵ) || throw(UndefKeywordError(:hϵ))  
    haskey(p, :tem∞) || throw(UndefKeywordError(:tem∞))  
    haskey(p, :energy) || throw(UndefKeywordError(:energy)) 

    kc::Float64 = float(p[:kc])
    kc ≤ 0 && throw(ArgumentError("`kc` should be larger than zero!\n"))
    aph::Float64 = float(p[:aph])
    aph ≤ 0 && throw(ArgumentError("`aph` should be larger than zero!\n"))
    cv::Float64 = float(p[:cv])
    cv ≤ 0 && throw(ArgumentError("`cv` should be larger than zero!\n"))
    rft::Float64 = float(p[:rft])
    h::Float64 = float(p[:h])
    h ≤ 0 && throw(ArgumentError("`h` should be larger than zero!\n"))
    hσ::Float64 = float(p[:hσ])
    hσ ≤ 0 && throw(ArgumentError("`hσ` should be larger than zero!\n"))
    hϵ::Float64 = float(p[:hϵ])
    hϵ ≤ 0 && throw(ArgumentError("`hϵ` should be larger than zero!\n"))
    tem∞::Float64 = float(p[:tem∞])
    energy::Float64 = float(p[:energy])
        
    return kc, aph, cv, rft, h, hσ, hϵ, tem∞, energy
end

@Peridynamics.params BBTTAMaterial BBTTAPointParameters

@Peridynamics.storage BBTTAMaterial struct BBTTAStorage <: Peridynamics.AbstractStorage
    @lthfield position::Matrix{Float64}
    @pointfield displacement::Matrix{Float64}
    @pointfield velocity::Matrix{Float64}
    @pointfield velocity_half::Matrix{Float64}
    @pointfield velocity_half_old::Matrix{Float64}
    @pointfield acceleration::Matrix{Float64}
    @pointfield b_int::Matrix{Float64}
    @pointfield b_int_old::Matrix{Float64}
    @pointfield b_ext::Matrix{Float64}
    @pointfield density_matrix::Matrix{Float64}
    @pointfield damage::Vector{Float64}
    bond_stretch::Vector{Float64}
    bond_active::Vector{Bool}
    @pointfield n_active_bonds::Vector{Int}
    @lthfield temperature::Matrix{Float64}
    @pointfield pflux::Matrix{Float64}
    @pointfield hsource::Matrix{Float64}
    @pointfield ablation_deep::Matrix{Float64}
    @lthfield exist::Matrix{Int} # 1: exist, 0: ablated
end

function Peridynamics.init_field_solver(::Thermstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:position})
    return copy(system.position)
end

function Peridynamics.init_field_solver(::Thermstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:displacement})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:velocity})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:velocity_half})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:acceleration})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:b_int})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:b_ext})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:temperature})
    return zeros(1, size(system.position, 2))
end

function Peridynamics.init_field_solver(::Thermstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:pflux})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:hsource})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:ablation_deep})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:exist})
    return ones(Int, 1, size(system.position, 2))
end

function Peridynamics.Peridynamics.init_field(::BBTTAMaterial, ::Peridynamics.AbstractTimeSolver, system::Peridynamics.BondSystem,
                    ::Val{:bond_stretch})
    return zeros(Peridynamics.get_n_bonds(system))
end

function Peridynamics.Peridynamics.allowed_material_kwargs(::BBTTAMaterial)
    return (thermal_temp_ablation_kwargs())
end

function calc_pflux_ablation!(chunk::Peridynamics.AbstractBodyChunk)
    (; system, mat, paramsetup, storage) = chunk
    storage.pflux .= 0.0
    storage.n_active_bonds .= 0
    for point_id in eachindex(chunk.system.chunk_handler.loc_points)
        pflux_point!(storage, system, mat, paramsetup, point_id)
    end
    return nothing
end

function pflux_point!(storage::BBTTAStorage, system::Peridynamics.BondSystem, 
    ::BBTTAMaterial, param::BBTTAPointParameters, i::Int) 
    
    if storage.exist[1, i] == 0
        storage.pflux[1, i] = 0.0
    else
        for bond_id in system.bond_ids[i]
            bond = system.bonds[bond_id]
            j, L = bond.neighbor, bond.length            
            Δtem = storage.temperature[1, j] - storage.temperature[1, i]
            kp_tm = 0.5*(param.kp_T(param.kc, param.kp, storage.temperature[1, i]) + param.kp_T(param.kc, param.kp, storage.temperature[1, j]))           
            storage.pflux[1, i] += storage.bond_active[bond_id] *  Δtem / L * kp_tm * system.volume[j] * storage.exist[1, j]
        end   
    end

    return nothing
end

function pflux_point!(storage::BBTTAStorage, system::Peridynamics.BondSystem, 
    ::BBTTAMaterial, paramhandler::Peridynamics.ParameterHandler, i::Int) 

    if storage.exist[1, i] == 0
        storage.pflux[1, i] = 0.0
    else
        params_i = Peridynamics.get_params(paramhandler, i)
        for bond_id in system.bond_ids[i]
            bond = system.bonds[bond_id]
            j, L = bond.neighbor, bond.length       
            Δtem = storage.temperature[1, j] - storage.temperature[1, i]
            params_j = Peridynamics.get_params(paramhandler, j)
            kp_tm = 0.5*(params_j.kp_T(params_j.kc, params_j.kp, storage.temperature[1, j]) + params_i.kp_T(params_i.kc, params_i.kp, storage.temperature[1, i]))s        
            storage.pflux[1, i] += storage.bond_active[bond_id] * kp_tm * Δtem / L * system.volume[j] * storage.exist[1, j]
        end
    end

    return nothing
end






