struct BBTAMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
    dmgmodel::DM
    function BBTAMaterial{C}(dmgmodel::DM) where {C,DM}
        new{C,DM}(dmgmodel)
    end
end

function BBTAMaterial{C}(; dmgmodel::Peridynamics.AbstractDamageModel=CriticalStretch()) where {C}
    return BBTAMaterial{C}(dmgmodel)
end
BBTAMaterial(; kwargs...) = BBTAMaterial{NoCorrection}(; kwargs...)

struct BBTAPointParameters <: Peridynamics.AbstractPointParameters
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

function BBTAPointParameters(mat::BBTAMaterial, p::Dict{Symbol,Any})
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

    if haskey(p, :thick) #1D_axis_symmetric 
        thick = float(p[:thick])
        bc = 9 * E / (π * thick * δ^3) # bond constant
        kp = kc / δ # microcndicitvity constant
    else #2D_axis_symmetric
        bc = 12 * E / (π * δ^4) 
        kp = 4 * kc / (π * δ^2) 
    end   
    return BBTAPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc, kc, kp, aph, cv, rft, h, hσ, hϵ, tem∞)
end


@Peridynamics.params BBTAMaterial BBTAPointParameters

@Peridynamics.storage BBTAMaterial struct BBTAStorage <: Peridynamics.AbstractStorage
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

function Peridynamics.Peridynamics.init_field(::BBTAMaterial, ::Peridynamics.AbstractTimeSolver, system::Peridynamics.BondSystem,
                    ::Val{:bond_stretch})
    return zeros(Peridynamics.get_n_bonds(system))
end

function Peridynamics.Peridynamics.allowed_material_kwargs(::BBTAMaterial)
    return (thermal_kwargs())
end

function pflux_point!(storage::BBTAStorage, system::Peridynamics.BondSystem, 
    ::BBTAMaterial, param::BBTAPointParameters, i::Int, mbd_t::Vector{Float64}) 

    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        mof_th = mbd_t[bond_id]
        vol_A = L * (param.δ/3) * π * (storage.position[1, j]+storage.position[1, i])
        
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]

        shape = 2*storage.position[1, j]/(storage.position[1, j] + storage.position[1, i]) 
        
        storage.pflux[1, i] += (param.kp * Δtem / L^2) * shape * vol_A  
    end
    return nothing
end

function pflux_point!(storage::BBTAStorage, system::Peridynamics.BondSystem, 
    ::BBTAMaterial, paramhandler::Peridynamics.ParameterHandler, i::Int, mbd_t::Vector{Float64}) 

    params_i = Peridynamics.get_params(paramhandler, i)
    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        mof_th = mbd_t[bond_id]

        shape = 2*storage.position[1, j]/(storage.position[1, j] + storage.position[1, i]) 
        
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]

        params_j = Peridynamics.get_params(paramhandler, j)
        
        storage.pflux[1, i] += storage.bond_active[bond_id] * (0.5*(params_i.kp + params_j.kp) * Δtem / L^2) * 
        system.volume[j] * shape
    end
    return nothing
end
