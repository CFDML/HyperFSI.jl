"""
Anisotropic module for thermal diffusion with temperature-dependent properties.
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

@inline thermal_temp_anis_kwargs() = (:E, :nu, :G, :K, :lambda, :mu, :rho, :horizon, 
:Gc, :epsilon_c, :kcx, :kcy, :aph, :cv, :rft, :h, :hσ, :hϵ, :tem∞, :thick, :kpx_T, :kpy_T, :cv_T, :aph_T)

struct BBTTxyMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
    dmgmodel::DM
    function BBTTxyMaterial{C}(dmgmodel::DM) where {C,DM} #BB_Temperature_dependent_Theermal
        new{C,DM}(dmgmodel)
    end
end

function BBTTxyMaterial{C}(; dmgmodel::Peridynamics.AbstractDamageModel=CriticalStretch()) where {C}
    return BBTTxyMaterial{C}(dmgmodel)
end

BBTTxyMaterial(; kwargs...) = BBTTxyMaterial{NoCorrection}(; kwargs...)

struct BBTTxyPointParameters <: Peridynamics.AbstractPointParameters
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
    kcx::Float64 # thermal conductivity at room temperature_x
    kcy::Float64 # thermal conductivity at room temperature_y  
    kpx::Float64 # microconductivity at room temperature_x
    kpy::Float64 # microconductivity at room temperature_y
    kp::Float64 # microconductivity at room temperature_maiximum   
    aph::Float64 # thermal expansion at room temperature
    cv::Float64 # specific heat capacity at room temperature
    rft::Float64 # Reference temperature at room temperature
    h::Float64 # convective heat transfer coefficient,
    hσ::Float64 # Stefan-Boltzman constant,
    hϵ::Float64 # emissivity
    tem∞::Float64 # temperature of the surrounding medium
    kp_Tx::Function # function describing kp variation with temperature f(kp,T)_x
    kp_Ty::Function # function describing kp variation with temperature f(kp,T)_y    
    cv_T::Function # function describing kp variation with temperature f(kp,T)
    aph_T::Function # function describing kp variation with temperature f(kp,T)
end

function BBTTxyPointParameters(mat::BBTTxyMaterial, p::Dict{Symbol,Any})
    par =  Peridynamics.get_given_elastic_params(p)
    (; E, nu, G, K, λ, μ) = par

    if haskey(p, :thick)
        p[:nu] = 1/3
    else 
        p[:nu] = 1/4    
    end
   
    (; δ, rho, E, nu, G, K, λ, μ) =  Peridynamics.get_required_point_parameters(mat, p)
    Gc, εc =  Peridynamics.get_frac_params(mat.dmgmodel, p, δ, K)
    kcx, kcy, aph, cv, rft, h, hσ, hϵ, tem∞ = get_thermal_anisotropic_params(p, δ)

    if haskey(p, :thick) #2D 
        thick = float(p[:thick])
        kpx = (6*kcy -2*kcx)/ (π *  δ^2) # microcndicitvity constant
        kpy = (8*kcx -8*kcy)/ (π *  δ^2) # microcndicitvity constant
        kp = max(kpx, kpy)  
        bc = 9 * E / (π * thick * δ^3)      
    else #3D
        error("Currently only 2D problems are supported!\n")
    end 

    function const_kpx(kcx, kpx, T)
        return kpx
    end

    function const_kpy(kcy, kpy, T)
        return kpy
    end

    function const_cv(cv, T)
        return cv
    end

    function const_aph(aph, T)
        return aph
    end   

    kp_Tx = get(p, :kp_Tx, const_kpx)
    kp_Ty = get(p, :kp_Ty, const_kpy)
    cv_T = get(p, :cv_T, const_cv)
    aph_T = get(p, :aph_T, const_aph)

    return BBTTxyPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc, kcx, kcy, kpx, kpy, kp, aph, cv, rft, h, hσ, hϵ, tem∞, kp_Tx, kp_Ty, cv_T, aph_T)
end

@Peridynamics.params BBTTxyMaterial BBTTxyPointParameters

@Peridynamics.storage BBTTxyMaterial struct BBTTxyStorage <: Peridynamics.AbstractStorage
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

function Peridynamics.Peridynamics.init_field(::BBTTxyMaterial, ::Peridynamics.AbstractTimeSolver, system::Peridynamics.BondSystem,
                    ::Val{:bond_stretch})
    return zeros(Peridynamics.get_n_bonds(system))
end

function Peridynamics.Peridynamics.allowed_material_kwargs(::BBTTxyMaterial)
    return (thermal_temp_anis_kwargs())
end

function pflux_point!(storage::BBTTxyStorage, system::Peridynamics.BondSystem, 
    ::BBTTxyMaterial, param::BBTTxyPointParameters, i::Int) 

    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]
        ae = ((storage.position[1, j] - storage.position[1, i]) / L)^2       

        kp_tmx = 0.5*(param.kp_Tx(param.kcx, param.kpx, storage.temperature[1, i]) + param.kp_Tx(param.kcx, param.kpx, storage.temperature[1, j]))
        kp_tmy = 0.5*(param.kp_Ty(param.kcy, param.kpy, storage.temperature[1, i]) + param.kp_Ty(param.kcy, param.kpy, storage.temperature[1, j]))  
        kp_tm = kp_tmx + kp_tmy * ae      
        storage.pflux[1, i] += storage.bond_active[bond_id] *  (Δtem / L^2) * kp_tm * system.volume[j]
    end
    return nothing
end

function pflux_point!(storage::BBTTxyStorage, system::Peridynamics.BondSystem, 
    ::BBTTxyMaterial, paramhandler::Peridynamics.ParameterHandler, i::Int) 

    params_i = Peridynamics.get_params(paramhandler, i)
    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length

        ae = ((storage.position[1, j] - storage.position[1, i]) / L)^2        
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]

        params_j = Peridynamics.get_params(paramhandler, j)

        kp_tmx = 0.5*(params_j.kp_Tx(params_j.kcx, params_j.kpx, storage.temperature[1, j]) + params_i.kp_Tx(params_i.kcx, params_i.kpx, storage.temperature[1, i]))
        kp_tmy = 0.5*(params_j.kp_Ty(params_j.kcy, params_j.kpy, storage.temperature[1, j]) + params_i.kp_Ty(params_i.kcy, params_i.kpy, storage.temperature[1, i]))  
        kp_tm = kp_tmx + kp_tmy * ae
        
        storage.pflux[1, i] += storage.bond_active[bond_id] * kp_tm * (Δtem / L^2) * system.volume[j]
    end
    return nothing
end

