"""
extends from Selda's work !! Not based on existing module
Anisotropic module for thermal diffusion with temperature-dependent properties.
State-Based Thermal Diffusion Module with Temperature-Dependent Properties

This module extends the standard bond-based thermal diffusion model
by incorporating temperature-dependent thermal properties, including:

Thermal conductivity tensor in classical theory (kc) 
Specific heat capacity (cv)
Thermal expansion coefficient (aph).

Customizable Functions:
kc_T(kc, T): 2-parameter function (classical Thermal conductivity, temperature T).
cv_T(cv, T): 2-parameter function (specific heat cv, temperature T).
aph_T(aph, T): 2-parameter function (expansion coefficient aph, temperature T).

Default Behavior:
If no function is assigned (e.g., kc_T is undefined), the property is treated as temperature-independent.
"""

@inline thermal_temp_anis_state_kwargs() = (:E, :nu, :G, :K, :lambda, :mu, :rho, :horizon, 
:Gc, :epsilon_c, :kc, :aph, :cv, :rft, :h, :hσ, :hϵ, :tem∞, :thick, :kc_T, :cv_T, :aph_T)

struct OSBTTMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
    dmgmodel::DM
    function OSBTTMaterial{C}(dmgmodel::DM) where {C,DM} # Ordinary_SB_Temperature_dependent_Theermal
        new{C,DM}(dmgmodel)
    end
end

function OSBTTMaterial{C}(; dmgmodel::Peridynamics.AbstractDamageModel=CriticalStretch()) where {C}
    return OSBTTMaterial{C}(dmgmodel)
end

OSBTTMaterial(; kwargs...) = OSBTTMaterial{NoCorrection}(; kwargs...)

struct OSBTTPointParameters <: Peridynamics.AbstractPointParameters
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
    kc::Matrix{Float64} # thermal conductivity tensor at room temperature 
    kp::Float64 # microconductivity at room temperature  
    ks::Float64 # element of shape tensor (approximation for non-uniform discretization)
    aph::Float64 # thermal expansion at room temperature
    cv::Float64 # specific heat capacity at room temperature
    rft::Float64 # Reference temperature at room temperature
    h::Float64 # convective heat transfer coefficient,
    hσ::Float64 # Stefan-Boltzman constant,
    hϵ::Float64 # emissivity
    tem∞::Float64 # temperature of the surrounding medium
    kc_T::Function # function describing kc variation with temperature f(kc,T) 
    cv_T::Function # function describing kp variation with temperature f(kp,T)
    aph_T::Function # function describing kp variation with temperature f(kp,T)
end

function OSBTTPointParameters(mat::OSBTTMaterial, p::Dict{Symbol,Any})
    par =  Peridynamics.get_given_elastic_params(p)
    (; E, nu, G, K, λ, μ) = par

    if haskey(p, :thick)
        p[:nu] = 1/3
    else 
        p[:nu] = 1/4    
    end
   
    (; δ, rho, E, nu, G, K, λ, μ) =  Peridynamics.get_required_point_parameters(mat, p)
    Gc, εc =  Peridynamics.get_frac_params(mat.dmgmodel, p, δ, K)
    kc, aph, cv, rft, h, hσ, hϵ, tem∞ = get_thermal_anisotropic_state_params(p, δ)
    
    kcm = maximum(kc) # using the maximum value of tensor for microconductivity
    if haskey(p, :thick) #2D 
        thick = float(p[:thick])
        bc = 9 * E / (π * thick * δ^3)   
        kp = (6*kcm) / (π * thick * δ^3) # microcndicitvity constant 
        ks = (2π * δ^4 * thick) / 3  
    else #3D
        bc = 12 * E / (π * δ^4)
        kp = (6*kcm) / (π * δ^4) # microcndicitvity constant
        ks = π * δ^5
    end 

    function const_kc(kc, T)
        return kc
    end

    function const_cv(cv, T)
        return cv
    end

    function const_aph(aph, T)
        return aph
    end   

    kc_T = get(p, :kc_T, const_kc)
    cv_T = get(p, :cv_T, const_cv)
    aph_T = get(p, :aph_T, const_aph)

    return OSBTTPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc, kc, kp, ks, aph, cv, rft, h, hσ, hϵ, tem∞, kc_T, cv_T, aph_T)
end

function get_thermal_anisotropic_state_params(p::Dict{Symbol,Any}, δ)
    haskey(p, :kc) || throw(UndefKeywordError(:kc))
    haskey(p, :aph) || throw(UndefKeywordError(:aph))
    haskey(p, :cv) || throw(UndefKeywordError(:cv))
    haskey(p, :rft) || throw(UndefKeywordError(:rft))
    haskey(p, :h) || throw(UndefKeywordError(:h))  
    haskey(p, :hσ) || throw(UndefKeywordError(:hσ))  
    haskey(p, :hϵ) || throw(UndefKeywordError(:hϵ))  
    haskey(p, :tem∞) || throw(UndefKeywordError(:tem∞))  

    kc::Matrix{Float64} = float.(p[:kc])
    if size(kc) != (3, 3)
        throw(ArgumentError("`kc` must be a 3×3 matrix, got $(size(kc))"))
    end
    all(kc .>= 0) || throw(ArgumentError("All elements in `kc` should be larger than zero!"))

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
        
    return kc, aph, cv, rft, h, hσ, hϵ, tem∞
end

@Peridynamics.params OSBTTMaterial OSBTTPointParameters

@Peridynamics.storage OSBTTMaterial struct OSBTTStorage <: Peridynamics.AbstractStorage
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
    @lthfield qc::Matrix{Float64} 
    @pointfield pflux::Matrix{Float64}
    @pointfield hsource::Matrix{Float64}
end


function Peridynamics.init_field_solver(::Thermstep_osb, system::Peridynamics.AbstractSystem, ::Val{:position})
    return copy(system.position)
end

function Peridynamics.init_field_solver(::Thermstep_osb, system::Peridynamics.AbstractSystem, ::Val{:displacement})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_osb, system::Peridynamics.AbstractSystem, ::Val{:velocity})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_osb, system::Peridynamics.AbstractSystem, ::Val{:velocity_half})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_osb, system::Peridynamics.AbstractSystem, ::Val{:acceleration})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_osb, system::Peridynamics.AbstractSystem, ::Val{:b_int})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_osb, system::Peridynamics.AbstractSystem, ::Val{:b_ext})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_osb, system::Peridynamics.AbstractSystem, ::Val{:temperature})
    return zeros(1, Peridynamics.size(system.position, 2))
end

function Peridynamics.init_field_solver(::Thermstep_osb, system::Peridynamics.AbstractSystem, ::Val{:qc})
    return zeros(3, Peridynamics.size(system.position, 2))
end

function Peridynamics.init_field_solver(::Thermstep_osb, system::Peridynamics.AbstractSystem, ::Val{:pflux})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_osb, system::Peridynamics.AbstractSystem, ::Val{:hsource})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.Peridynamics.init_field(::OSBTTMaterial, ::Peridynamics.AbstractTimeSolver, system::Peridynamics.BondSystem,
                    ::Val{:bond_stretch})
    return zeros(Peridynamics.get_n_bonds(system))
end

function Peridynamics.Peridynamics.allowed_material_kwargs(::OSBTTMaterial)
    return (thermal_temp_anis_state_kwargs())
end

function qc_density_point!(storage::OSBTTStorage, system::Peridynamics.BondSystem, 
    ::OSBTTMaterial, param::OSBTTPointParameters, t, Δt, i::Int) 

    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]
        ae = (storage.position[:, j] - storage.position[:, i]) ./ L 

        kc_tm = 0.5.*(param.kc_T(param.kc, storage.temperature[1, i]) + 
                      param.kc_T(param.kc, storage.temperature[1, j]))             
    
        storage.qc[:, i] += kc_tm * (Δtem * system.volume[j] .* ae) * param.δ / param.ks
    end       
    return nothing
end

function qc_density_point!(storage::OSBTTStorage, system::Peridynamics.BondSystem, 
    ::OSBTTMaterial, paramhandler::Peridynamics.ParameterHandler, t, Δt, i::Int) 

    param_i = Peridynamics.get_params(paramhandler, i)
    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length

        Δtem = storage.temperature[1, j] - storage.temperature[1, i]
        ae = (storage.position[:, j] - storage.position[:, i]) ./ L      

        param_j = Peridynamics.get_params(paramhandler, j)
        kc_tm = 0.5.*(param_i.kc_T(param_i.kc, storage.temperature[1, i]) + 
                      param_j.kc_T(param_j.kc, storage.temperature[1, j])) 
        
        storage.qc[:, i] += kc_tm * (Δtem * system.volume[j] .* ae) * param.δ / param.ks
    end
    return nothing
end

function pflux_density_point!(storage::OSBTTStorage, system::Peridynamics.BondSystem, 
    ::OSBTTMaterial, param::OSBTTPointParameters, t, Δt, i::Int) 

    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
         
        qc = storage.qc[:, i] + storage.qc[:, j] 

        qflux = dot(qc, (storage.position[:, j] - storage.position[:, i])) / L * param.δ / param.ks
        storage.pflux[1, i] += storage.bond_active[bond_id] *  qflux * system.volume[j]
    end
    return nothing
end

function pflux_density_point!(storage::OSBTTStorage, system::Peridynamics.BondSystem, 
    ::OSBTTMaterial, paramhandler::Peridynamics.ParameterHandler, t, Δt, i::Int) 

    param_i = Peridynamics.get_params(paramhandler, i)
    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        qc = storage.qc[:, i] + storage.qc[:, j] 
    
        qflux = dot(qc, (storage.position[:, j] - storage.position[:, i])) / L * param_i.δ / param_i.ks
        storage.pflux[1, i] += storage.bond_active[bond_id] *  qflux * system.volume[j]
    end
    return nothing
end

function calc_qc_density!(dh::Peridynamics.ThreadsBodyDataHandler, t, Δt)
    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_loc_to_halo!(dh, chunk_id)
        calc_qc_density!(dh.chunks[chunk_id], t, Δt)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_halo_to_loc!(dh, chunk_id)
    end
    return nothing
end

function calc_qc_density!(chunk::Peridynamics.AbstractBodyChunk{<:Peridynamics.AbstractBondSystem}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    storage.qc .= 0
    for point_id in Peridynamics.each_point_idx(chunk)
        qc_density_point!(storage, system, mat, paramsetup, t, Δt, point_id)
    end
    nancheck_osb(chunk, t)
    return nothing
end

function calc_pflux_density!(dh::Peridynamics.ThreadsBodyDataHandler, t, Δt)
    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_loc_to_halo!(dh, chunk_id)
        calc_pflux_density!(dh.chunks[chunk_id], t, Δt)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_halo_to_loc!(dh, chunk_id)
    end
    return nothing
end

function calc_pflux_density!(chunk::Peridynamics.AbstractBodyChunk{<:Peridynamics.AbstractBondSystem}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    storage.pflux .= 0
    for point_id in Peridynamics.each_point_idx(chunk)
        pflux_density_point!(storage, system, mat, paramsetup, t, Δt, point_id)
    end
    nancheck_osb(chunk, t)
    return nothing
end

function nancheck_osb(chunk::Peridynamics.AbstractBodyChunk, t)
    if Peridynamics.containsnan(chunk.storage.qc) || Peridynamics.containsnan(chunk.storage.pflux)
        msg = "NaN's found in field `qc or pflux` at simulation time $(t)!\n"
        error(msg)
    end
    return nothing
end
