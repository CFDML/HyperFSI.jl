@inline thermal_kwargs() = (:E, :nu, :G, :K, :lambda, :mu, :rho, :horizon, 
:Gc, :epsilon_c, :kc, :aph, :cv, :rft, :h, :hσ, :hϵ, :tem∞, :thick)

struct BBTMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
    dmgmodel::DM
    function BBTMaterial{C}(dmgmodel::DM) where {C,DM}
        new{C,DM}(dmgmodel)
    end
end

function BBTMaterial{C}(; dmgmodel::Peridynamics.AbstractDamageModel=CriticalStretch()) where {C}
    return BBTMaterial{C}(dmgmodel)
end
BBTMaterial(; kwargs...) = BBTMaterial{NoCorrection}(; kwargs...)

struct BBTPointParameters <: Peridynamics.AbstractPointParameters
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
    kc::Float64 # thermal conductivity
    kp::Float64 # microconductivity
    aph::Float64 # thermal expansion
    cv::Float64 # specific heat capacity
    rft::Float64 # Reference temperature
    h::Float64 # convective heat transfer coefficient,
    hσ::Float64 # Stefan-Boltzman constant,
    hϵ::Float64 # emissivity
    tem∞::Float64 # temperature of the surrounding medium
end

function BBTPointParameters(mat::BBTMaterial, p::Dict{Symbol,Any})
    par = Peridynamics.get_given_elastic_params(p)
    (; E, nu, G, K, λ, μ) = par

    if haskey(p, :thick)
        p[:nu] = 1/3
    else 
        p[:nu] = 1/4    
    end
   
    (; δ, rho, E, nu, G, K, λ, μ) = Peridynamics.get_required_point_parameters(mat, p)
    Gc, εc = Peridynamics.get_frac_params(mat.dmgmodel, p, δ, K)
    kc, aph, cv, rft, h, hσ, hϵ, tem∞ = get_thermal_params(p, δ)

    if haskey(p, :thick) #2D 
        thick = float(p[:thick])
        bc = 9 * E / (π * thick * δ^3) # bond constant
        kp = 6 * kc / (π * thick * δ^3) # microcndicitvity constant
    else #3D
        bc = 12 * E / (π * δ^4) 
        kp = 6 * kc / (π * δ^4) 
    end   
    return BBTPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc, kc, kp, aph, cv, rft, h, hσ, hϵ, tem∞)
end

function get_thermal_params(p::Dict{Symbol,Any}, δ)
    haskey(p, :kc) || throw(UndefKeywordError(:kc))
    haskey(p, :aph) || throw(UndefKeywordError(:aph))
    haskey(p, :cv) || throw(UndefKeywordError(:cv))
    haskey(p, :rft) || throw(UndefKeywordError(:rft))
    haskey(p, :h) || throw(UndefKeywordError(:h))  
    haskey(p, :hσ) || throw(UndefKeywordError(:hσ))  
    haskey(p, :hϵ) || throw(UndefKeywordError(:hϵ))  
    haskey(p, :tem∞) || throw(UndefKeywordError(:tem∞))  

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
        
    return kc, aph, cv, rft, h, hσ, hϵ, tem∞
end

@Peridynamics.params BBTMaterial BBTPointParameters

@Peridynamics.storage BBTMaterial struct BBTStorage <: Peridynamics.AbstractStorage
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

function Peridynamics.init_field_solver(::Thermstep, system::Peridynamics.AbstractSystem, ::Val{:position})
    return copy(system.position)
end

function Peridynamics.init_field_solver(::Thermstep, system::Peridynamics.AbstractSystem, ::Val{:displacement})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep, system::Peridynamics.AbstractSystem, ::Val{:velocity})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep, system::Peridynamics.AbstractSystem, ::Val{:velocity_half})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep, system::Peridynamics.AbstractSystem, ::Val{:acceleration})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep, system::Peridynamics.AbstractSystem, ::Val{:b_int})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep, system::Peridynamics.AbstractSystem, ::Val{:b_ext})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep, system::Peridynamics.AbstractSystem, ::Val{:temperature})
    return zeros(1, Peridynamics.size(system.position, 2))
end

function Peridynamics.init_field_solver(::Thermstep, system::Peridynamics.AbstractSystem, ::Val{:pflux})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermstep, system::Peridynamics.AbstractSystem, ::Val{:hsource})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.Peridynamics.init_field(::BBTMaterial, ::Peridynamics.AbstractTimeSolver, system::Peridynamics.BondSystem,
                    ::Val{:bond_stretch})
    return zeros(Peridynamics.get_n_bonds(system))
end

function Peridynamics.Peridynamics.allowed_material_kwargs(::BBTMaterial)
    return (thermal_kwargs())
end

# Thermal 
@inline function update_temperature!(b::Peridynamics.AbstractBodyChunk, Δt::Float64)
    for point_id in Peridynamics.each_point_idx(b)
        param = Peridynamics.get_params(b, point_id)
        k = Δt / (param.rho * param.cv)
        _update_temperature!(b.storage.temperature, b.storage.pflux, b.storage.hsource, k, point_id)
    end
    return nothing
end

@inline function _update_temperature!(temperature, pflux, hsource, k, i)
    temperature[1, i] += (pflux[1, i] + hsource[1, i]) * k
    return nothing
end

function calc_pflux!(chunk::Peridynamics.AbstractBodyChunk, mbd_t::Vector{Float64})
    (; system, mat, paramsetup, storage) = chunk
    storage.pflux .= 0.0
    storage.n_active_bonds .= 0
    for point_id in eachindex(chunk.system.chunk_handler.loc_points)
        pflux_point!(storage, system, mat, paramsetup, point_id, mbd_t)
    end
    return nothing
end

function pflux_point!(storage::BBTStorage, system::Peridynamics.BondSystem, 
    ::BBTMaterial, param::BBTPointParameters, i::Int, mbd_t::Vector{Float64}) 

    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        mof_th = mbd_t[bond_id]
        
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]
        
        storage.pflux[1, i] += storage.bond_active[bond_id] * (param.kp * Δtem / L) * 
        system.volume[j] #* mof_th
    end
    return nothing
end

function pflux_point!(storage::BBTStorage, system::Peridynamics.BondSystem, 
    ::BBTMaterial, paramhandler::Peridynamics.ParameterHandler, i::Int, mbd_t::Vector{Float64}) 

    params_i = Peridynamics.get_params(paramhandler, i)
    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        mof_th = mbd_t[bond_id]
        
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]

        params_j = Peridynamics.get_params(paramhandler, j)
        
        storage.pflux[1, i] += storage.bond_active[bond_id] * (0.5*(params_i.kp + params_j.kp) * Δtem / L) * 
        system.volume[j] * mof_th
    end
    return nothing
end
