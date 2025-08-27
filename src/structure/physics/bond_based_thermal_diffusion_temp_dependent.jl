"""
Bond-Based Thermal Diffusion Module with Temperature-Dependent Properties

This module extends the standard bond-based thermal diffusion model
by incorporating temperature-dependent thermal properties, including:

Thermal conductivity in PD (kp) 
Specific heat capacity (cv)
Thermal expansion coefficient (aph).

Customizable Functions:
kp_T(kc, kp, T): 3-parameter function (classical Thermal conductivity, PD Thermal conductivity , temperature T).
cv_T(cv, T): 2-parameter function (specific heat cv, temperature T).
aph_T(aph, T): 2-parameter function (expansion coefficient aph, temperature T).

Default Behavior:
If no function is assigned (e.g., kp_T is undefined), the property is treated as temperature-independent.
"""

@inline thermal_temp_kwargs() = (:E, :nu, :G, :K, :lambda, :mu, :rho, :horizon, 
:Gc, :epsilon_c, :kc, :aph, :cv, :rft, :h, :hσ, :hϵ, :tem∞, :thick, :kp_T, :cv_T, :aph_T)

struct BBTTMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
    dmgmodel::DM
    function BBTTMaterial{C}(dmgmodel::DM) where {C,DM} #BB_Temperature_dependent_Theermal
        new{C,DM}(dmgmodel)
    end
end

function BBTTMaterial{C}(; dmgmodel::Peridynamics.AbstractDamageModel=CriticalStretch()) where {C}
    return BBTTMaterial{C}(dmgmodel)
end

BBTTMaterial(; kwargs...) = BBTTMaterial{NoCorrection}(; kwargs...)

struct BBTTPointParameters <: Peridynamics.AbstractPointParameters
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
end

function const_kp(kc, kp, T)
    return kp
end

function const_cv(cv, T)
    return cv
end

function const_aph(aph, T)
    return aph
end

function BBTTPointParameters(mat::BBTTMaterial, p::Dict{Symbol,Any})
    par =  Peridynamics.get_given_elastic_params(p)
    (; E, nu, G, K, λ, μ) = par

    if haskey(p, :thick)
        p[:nu] = 1/3
    else 
        p[:nu] = 1/4    
    end
   
    (; δ, rho, E, nu, G, K, λ, μ) =  Peridynamics.get_required_point_parameters(mat, p)
    Gc, εc =  Peridynamics.get_frac_params(mat.dmgmodel, p, δ, K)
    kc, aph, cv, rft, h, hσ, hϵ, tem∞ = get_thermal_params(p, δ)

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

    return BBTTPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc, kc, kp, aph, cv, rft, h, hσ, hϵ, tem∞, kp_T, cv_T, aph_T)
end

@Peridynamics.params BBTTMaterial BBTTPointParameters

@Peridynamics.storage BBTTMaterial struct BBTTStorage <: Peridynamics.AbstractStorage
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
end

function Peridynamics.Peridynamics.init_field(::BBTTMaterial, ::Peridynamics.AbstractTimeSolver, system::Peridynamics.BondSystem,
                    ::Val{:bond_stretch})
    return zeros(Peridynamics.get_n_bonds(system))
end

function Peridynamics.Peridynamics.allowed_material_kwargs(::BBTTMaterial)
    return (thermal_temp_kwargs())
end

function pflux_point!(storage::BBTTStorage, system::Peridynamics.BondSystem, 
    ::BBTTMaterial, param::BBTTPointParameters, i::Int, mbd_t::Vector{Float64}) 

    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        mof_th = mbd_t[bond_id]
        
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]

        kp_tm = 0.5*(param.kp_T(param.kc, param.kp, storage.temperature[1, i]) + param.kp_T(param.kc, param.kp, storage.temperature[1, j]))
        
        storage.pflux[1, i] += storage.bond_active[bond_id] *  Δtem / L * kp_tm * system.volume[j] * mof_th
    end
    return nothing
end

function pflux_point!(storage::BBTTStorage, system::Peridynamics.BondSystem, 
    ::BBTTMaterial, paramhandler::Peridynamics.ParameterHandler, i::Int, mbd_t::Vector{Float64}) 

    params_i = Peridynamics.get_params(paramhandler, i)
    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        mof_th = mbd_t[bond_id]
        
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]

        params_j = Peridynamics.get_params(paramhandler, j)

        kp_tm = 0.5*(params_j.kp_T(params_j.kc, params_j.kp, storage.temperature[1, j]) + params_i.kp_T(params_i.kc, params_i.kp, storage.temperature[1, i]))
        
        storage.pflux[1, i] += storage.bond_active[bond_id] * kp_tm * Δtem / L * system.volume[j] * mof_th
    end
    return nothing
end
