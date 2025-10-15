"""
Anisotropic bond-based peridynamic model for thermo-mechanical problems.
Different thermal expansion coefficients in different directions can be specified.
Different with Selda's work, t_bond is applied.
Currently only support 2D problems.
"""
@inline thermal_anisotropic_kwargs() = (:E, :nu, :G, :K, :lambda, :mu, :rho, :horizon, 
:Gc, :epsilon_c, :kcx, :kcy, :aph, :cv, :rft, :h, :hσ, :hϵ, :tem∞, :thick)

struct BBTxyMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
    dmgmodel::DM
    function BBTxyMaterial{C}(dmgmodel::DM) where {C,DM}
        new{C,DM}(dmgmodel)
    end
end

function BBTxyMaterial{C}(; dmgmodel::Peridynamics.AbstractDamageModel=CriticalStretch()) where {C}
    return BBTxyMaterial{C}(dmgmodel)
end
BBTxyMaterial(; kwargs...) = BBTxyMaterial{NoCorrection}(; kwargs...)

struct BBTxyPointParameters <: Peridynamics.AbstractPointParameters
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
    kcx::Float64 # thermal conductivity
    kcy::Float64 # thermal conductivity
    kpx::Float64 # microconductivity
    kpy::Float64 # microconductivity   
    kp::Float64 # microconductivity 
    aph::Float64 # thermal expansion
    cv::Float64 # specific heat capacity
    rft::Float64 # Reference temperature
    h::Float64 # convective heat transfer coefficient,
    hσ::Float64 # Stefan-Boltzman constant,
    hϵ::Float64 # emissivity
    tem∞::Float64 # temperature of the surrounding medium
end

function BBTxyPointParameters(mat::BBTxyMaterial, p::Dict{Symbol,Any})
    par = Peridynamics.get_given_elastic_params(p)
    (; E, nu, G, K, λ, μ) = par

    if haskey(p, :thick)
        p[:nu] = 1/3
    else 
        p[:nu] = 1/4    
    end
   
    (; δ, rho, E, nu, G, K, λ, μ) = Peridynamics.get_required_point_parameters(mat, p)
    Gc, εc = Peridynamics.get_frac_params(mat.dmgmodel, p, δ, K)
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
    return BBTxyPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc, kcx, kcy, kpx, kpy, kp, aph, cv, rft, h, hσ, hϵ, tem∞)
end

function get_thermal_anisotropic_params(p::Dict{Symbol,Any}, δ)
    haskey(p, :kcx) || throw(UndefKeywordError(:kcx))
    haskey(p, :kcy) || throw(UndefKeywordError(:kcy))
    haskey(p, :aph) || throw(UndefKeywordError(:aph))
    haskey(p, :cv) || throw(UndefKeywordError(:cv))
    haskey(p, :rft) || throw(UndefKeywordError(:rft))
    haskey(p, :h) || throw(UndefKeywordError(:h))  
    haskey(p, :hσ) || throw(UndefKeywordError(:hσ))  
    haskey(p, :hϵ) || throw(UndefKeywordError(:hϵ))  
    haskey(p, :tem∞) || throw(UndefKeywordError(:tem∞))  

    kcx::Float64 = float(p[:kcx])
    kcx ≤ 0 && throw(ArgumentError("`kcx` should be larger than zero!\n"))
    kcy::Float64 = float(p[:kcy])
    kcy ≤ 0 && throw(ArgumentError("`kcy` should be larger than zero!\n"))
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
        
    return kcx, kcy, aph, cv, rft, h, hσ, hϵ, tem∞
end

@Peridynamics.params BBTxyMaterial BBTxyPointParameters

@Peridynamics.storage BBTxyMaterial struct BBTxyStorage <: Peridynamics.AbstractStorage
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

function Peridynamics.init_field_solver(::Thermstep_anis, system::Peridynamics.AbstractSystem, ::Val{:position})
    return copy(system.position)
end

function Peridynamics.init_field_solver(::Thermstep_anis, system::Peridynamics.AbstractSystem, ::Val{:displacement})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_anis, system::Peridynamics.AbstractSystem, ::Val{:velocity})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_anis, system::Peridynamics.AbstractSystem, ::Val{:velocity_half})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_anis, system::Peridynamics.AbstractSystem, ::Val{:acceleration})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_anis, system::Peridynamics.AbstractSystem, ::Val{:b_int})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_anis, system::Peridynamics.AbstractSystem, ::Val{:b_ext})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_anis, system::Peridynamics.AbstractSystem, ::Val{:temperature})
    return zeros(1, Peridynamics.size(system.position, 2))
end

function Peridynamics.init_field_solver(::Thermstep_anis, system::Peridynamics.AbstractSystem, ::Val{:pflux})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep_anis, system::Peridynamics.AbstractSystem, ::Val{:hsource})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.Peridynamics.init_field(::BBTxyMaterial, ::Peridynamics.AbstractTimeSolver, system::Peridynamics.BondSystem,
                    ::Val{:bond_stretch})
    return zeros(Peridynamics.get_n_bonds(system))
end

function Peridynamics.Peridynamics.allowed_material_kwargs(::BBTxyMaterial)
    return (thermal_anisotropic_kwargs())
end


function pflux_point!(storage::BBTxyStorage, system::Peridynamics.BondSystem, 
    ::BBTxyMaterial, param::BBTxyPointParameters, i::Int) 

    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]
        ae = ((storage.position[1, j] - storage.position[1, i]) / L)^2
        
        storage.pflux[1, i] += storage.bond_active[bond_id] * ((param.kpx + param.kpy * ae) * Δtem / L^2) * 
        system.volume[j]
    end
    return nothing
end

function pflux_point!(storage::BBTxyStorage, system::Peridynamics.BondSystem, 
    ::BBTxyMaterial, paramhandler::Peridynamics.ParameterHandler, i::Int) 

    params_i = Peridynamics.get_params(paramhandler, i)
    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        ae = ((storage.position[1, j] - storage.position[1, i]) / L)^2
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]

        params_j = Peridynamics.get_params(paramhandler, j)
        
        k = (params_i.kpx + params_j.kpx + (params_i.kpy + params_j.kpy) * ae) / 2
        storage.pflux[1, i] += storage.bond_active[bond_id] * (k * Δtem / L^2) * system.volume[j]
    end
    return nothing
end
