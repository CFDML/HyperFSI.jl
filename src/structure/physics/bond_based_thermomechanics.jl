@inline thermomech_kwargs() = (:E, :nu, :G, :K, :lambda, :mu, :rho, :horizon, 
:Gc, :epsilon_c, :kc, :aph, :cv, :rft, :h, :hσ, :hϵ, :tem∞, :thick)

struct BBTMMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
    dmgmodel::DM
    function BBTMMaterial{C}(dmgmodel::DM) where {C,DM}
        new{C,DM}(dmgmodel)
    end
end

function BBTMMaterial{C}(; dmgmodel::Peridynamics.AbstractDamageModel=CriticalStretch()) where {C}
    return BBTMMaterial{C}(dmgmodel)
end
BBTMMaterial(; kwargs...) = BBTMMaterial{NoCorrection}(; kwargs...)

struct BBTMPointParameters <: Peridynamics.AbstractPointParameters
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

function BBTMPointParameters(mat::BBTMMaterial, p::Dict{Symbol,Any})
    par = Peridynamics.get_given_elastic_params(p)
    (; E, nu, G, K, λ, μ) = par

    if haskey(p, :thick)
        p[:nu] = 1/3
    else 
        p[:nu] = 1/4    
    end
   
    (; δ, rho, E, nu, G, K, λ, μ) = Peridynamics.get_required_point_parameters(mat, p)
    Gc, εc = Peridynamics.get_frac_params(mat.dmgmodel, p, δ, K)
    kc, aph, cv, rft, h, hσ, hϵ, tem∞ = get_thermomech_params(p, δ)

    if haskey(p, :thick) #2D 
        thick = float(p[:thick])
        bc = 9 * E / (π * thick * δ^3) # bond constant
        kp = 6 * kc / (π * thick * δ^3) # microcndicitvity constant
    else #3D
        bc = 12 * E / (π * δ^4) 
        kp = 6 * kc / (π * δ^4) 
    end   
    return BBTMPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc, kc, kp, aph, cv, rft, h, hσ, hϵ, tem∞)
end

function get_thermomech_params(p::Dict{Symbol,Any}, δ)
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

@Peridynamics.params BBTMMaterial BBTMPointParameters

@Peridynamics.storage BBTMMaterial struct BBTMStorage <: Peridynamics.AbstractStorage
    @lthfield position::Matrix{Float64}
    @pointfield displacement::Matrix{Float64}
    @lthfield velocity::Matrix{Float64}
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


function Peridynamics.init_field_solver(::Thermomechstep, system::Peridynamics.AbstractSystem, ::Val{:position})
    return copy(system.position)
end

function Peridynamics.init_field_solver(::Thermomechstep, system::Peridynamics.AbstractSystem, ::Val{:displacement})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep, system::Peridynamics.AbstractSystem, ::Val{:velocity})
    return zeros(3, size(system.position, 2))
end

function Peridynamics.init_field_solver(::Thermomechstep, system::Peridynamics.AbstractSystem, ::Val{:velocity_half})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep, system::Peridynamics.AbstractSystem, ::Val{:acceleration})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep, system::Peridynamics.AbstractSystem, ::Val{:b_int})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep, system::Peridynamics.AbstractSystem, ::Val{:b_ext})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep, system::Peridynamics.AbstractSystem, ::Val{:temperature})
    return zeros(1, size(system.position, 2))
end

function Peridynamics.init_field_solver(::Thermomechstep, system::Peridynamics.AbstractSystem, ::Val{:pflux})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep, system::Peridynamics.AbstractSystem, ::Val{:hsource})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.Peridynamics.init_field(::BBTMMaterial, ::Peridynamics.AbstractTimeSolver, system::Peridynamics.BondSystem,
                    ::Val{:bond_stretch})
    return zeros(Peridynamics.get_n_bonds(system))
end

function Peridynamics.Peridynamics.allowed_material_kwargs(::BBTMMaterial)
    return (thermomech_kwargs())
end

function new_calc_force_density!(dh::Peridynamics.ThreadsBodyDataHandler, m_factors::Vector{Vector{Float64}}, t, Δt)
    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_loc_to_halo!(dh, chunk_id)
        new_calc_force_density!(dh.chunks[chunk_id], m_factors[chunk_id], t, Δt)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_halo_to_loc!(dh, chunk_id)
    end
    return nothing
end

function new_calc_force_density!(chunk::Peridynamics.AbstractBodyChunk{<:Peridynamics.AbstractBondSystem}, mbd_m::Vector{Float64}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    (; dmgmodel) = mat
    storage.b_int .= 0
    storage.n_active_bonds .= 0
    for point_id in Peridynamics.each_point_idx(chunk)
        cou_calc_failure!(storage, system, mat, dmgmodel, paramsetup, point_id)
        Peridynamics.calc_damage!(storage, system, mat, dmgmodel, paramsetup, point_id)
        new_force_density_point!(storage, system, mat, paramsetup, point_id, mbd_m)
    end
    Peridynamics.nancheck(chunk, t, Δt)
    return nothing
end

function cou_calc_failure!(storage::BBTMStorage, system::Peridynamics.BondSystem,
                       ::BBTMMaterial, ::Peridynamics.CriticalStretch,
                       paramsetup::BBTMPointParameters, i)
    (; εc, aph) = Peridynamics.get_params(paramsetup, i)
    (; position, n_active_bonds, bond_active, bond_stretch) = storage
    (; bonds) = system
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = Peridynamics.get_vector_diff(position, i, j)
        l = norm(Δxij)
        ε = (l - L) / L - 0.5 * (storage.temperature[1, j] + storage.temperature[1, i]) * aph
        bond_stretch[bond_id] = ε / l # note that this is  ε / l!
        if ε > εc && bond.fail_permit
            bond_active[bond_id] = false
        end
        n_active_bonds[i] += bond_active[bond_id]
    end
    return nothing
end

function cou_calc_failure!(storage::BBTMStorage, system::Peridynamics.BondSystem,
                       ::BBTMMaterial, ::Peridynamics.CriticalStretch,
                       paramhandler::Peridynamics.ParameterHandler, i)
    params_i = Peridynamics.get_params(paramhandler, i)
    (; position, n_active_bonds, bond_active, bond_stretch) = storage
    (; bonds) = system
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        params_j = Peridynamics.get_params(paramhandler, j)
        Δxij = Peridynamics.get_vector_diff(position, i, j)
        l = norm(Δxij)
        ε = (l - L) / L - 0.5 * (storage.temperature[1, j] + storage.temperature[1, i]) * (params_i.aph + params_j.aph) / 2
        bond_stretch[bond_id] = ε / l # note that this is  ε / l!
        εcm = min(params_i.εc, params_j.εc)^2 / max(params_i.εc, params_j.εc)
        if ε > εcm && bond.fail_permit
            bond_active[bond_id] = false
        end
        n_active_bonds[i] += bond_active[bond_id]
    end
    return nothing
end

function new_force_density_point!(storage::BBTMStorage, system::Peridynamics.BondSystem, ::BBTMMaterial,
                              params::BBTMPointParameters, i, mbd_m::Vector{Float64})
    (; position, bond_stretch, bond_active, b_int) = storage
    (; bonds, volume) = system
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        Δxij = Peridynamics.get_vector_diff(position, i, j)
        ε = bond_stretch[bond_id]
        ω = bond_active[bond_id] #* 1.0 * mbd_m[bond_id]
        b = ω * params.bc * ε * volume[j] .* Δxij
        Peridynamics.update_add_vector!(b_int, i, b)
    end
    return nothing
end

function new_force_density_point!(storage::BBTMStorage, system::Peridynamics.BondSystem, ::BBTMMaterial,
                              paramhandler::Peridynamics.ParameterHandler, i, mbd_m::Vector{Float64})
    (; position, bond_stretch, bond_active, b_int) = storage
    (; bonds, volume) = system
    params_i = Peridynamics.get_params(paramhandler, i)
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        Δxij = Peridynamics.get_vector_diff(position, i, j)
        ε = bond_stretch[bond_id]
        params_j = Peridynamics.get_params(paramhandler, j)
        ω = bond_active[bond_id] #* mbd_m[bond_id] 
        b = ω * (params_i.bc + params_j.bc) / 2 * ε * volume[j] .* Δxij
        Peridynamics.update_add_vector!(b_int, i, b)
    end
    return nothing
end

function cou_calc_pflux!(chunk::Peridynamics.AbstractBodyChunk, mbd_t::Vector{Float64})
    (; system, mat, paramsetup, storage) = chunk
    storage.pflux .= 0.0
    for point_id in eachindex(chunk.system.chunk_handler.loc_points)
        cou_pflux_point!(storage, system, mat, paramsetup, point_id, mbd_t)
    end
    return nothing
end

function cou_pflux_point!(storage::BBTMStorage, system::Peridynamics.BondSystem, 
    ::BBTMMaterial, param::BBTMPointParameters, i::Int, mbd_t::Vector{Float64}) 

    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        mof_th = mbd_t[bond_id]

        Δxijx = storage.position[1, j] - storage.position[1, i]
        Δxijy = storage.position[2, j] - storage.position[2, i]
        Δxijz = storage.position[3, j] - storage.position[3, i]
        Δvijx = storage.velocity[1, j] - storage.velocity[1, i]
        Δvijy = storage.velocity[2, j] - storage.velocity[2, i]
        Δvijz = storage.velocity[3, j] - storage.velocity[3, i]

        l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)
        ev = (Δxijx * Δvijx + Δxijy * Δvijy + Δxijz * Δvijz) / l
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]
        storage.pflux[1, i] += storage.bond_active[bond_id] * (param.kp * Δtem / L -
                    ev * 0.5 * param.rft * param.bc * param.aph) * system.volume[j] #* mof_th
    end
    return nothing
end

function cou_pflux_point!(storage::BBTMStorage, system::Peridynamics.BondSystem, 
    ::BBTMMaterial, paramhandler::Peridynamics.ParameterHandler, i::Int, mbd_t::Vector{Float64}) 

    params_i = Peridynamics.get_params(paramhandler, i)
    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        mof_th = 1.0 #mbd_t[bond_id]

        Δxijx = storage.position[1, j] - storage.position[1, i]
        Δxijy = storage.position[2, j] - storage.position[2, i]
        Δxijz = storage.position[3, j] - storage.position[3, i]
        Δvijx = storage.velocity[1, j] - storage.velocity[1, i]
        Δvijy = storage.velocity[2, j] - storage.velocity[2, i]
        Δvijz = storage.velocity[3, j] - storage.velocity[3, i]

        l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)
        ev = (Δxijx * Δvijx + Δxijy * Δvijy + Δxijz * Δvijz) / l
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]

        params_j = Peridynamics.get_params(paramhandler, j)

        storage.pflux[1, i] += storage.bond_active[bond_id] * (0.5 * (params_i.kp + params_j.kp) * Δtem / L - ev * 0.5 * 0.5 * (params_i.rft + params_j.rft) 
                            * 0.5 * (params_i.bc + params_j.bc) *  0.5 * (params_i.aph + params_j.aph)) * system.volume[j] * mof_th
    end
    return nothing
end

function ccalc_force_density!(chunk::Peridynamics.AbstractBodyChunk, mbd_m::Vector{Float64}, t)
    (; system, mat, paramsetup, storage) = chunk
    storage.b_int .= 0
    storage.n_active_bonds .= 0
    for point_id in Peridynamics.each_point_idx(chunk)
        cou_force_density_point!(storage, system, mat, paramsetup, point_id, mbd_m)
    end
    Peridynamics.nancheck(chunk, t)
    return nothing
end

function cou_force_density_point!(storage::BBTMStorage, system::Peridynamics.BondSystem, 
    ::BBTMMaterial, params::BBTMPointParameters, i::Int, mbd_m::Vector{Float64}) 

    (; εc, aph, bc) = Peridynamics.get_params(params, i)
    (; position, n_active_bonds, bond_active, temperature, b_int) = storage
    (; bonds, volume) = system
    
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length 
        mof_mh =  1.0 #mbd_m[bond_id]

        Δxij = Peridynamics.get_vector_diff(position, i, j)
        temc_avg = 0.5 * (temperature[1, j] + temperature[1, i])

        l = sqrt(Δxij[1] * Δxij[1] + Δxij[2] * Δxij[2] + Δxij[3] * Δxij[3])

        εm = (l - L) / L
        εt = temc_avg * aph
        ε = εm - εt

        if ε > εc && bond.fail_permit
            bond_active[bond_id] = false
        end
        n_active_bonds[i] += bond_active[bond_id]
        
        b_n = bond_active[bond_id] * bc * ε/l * volume[j] .* Δxij * mof_mh
        b_int[:, i] += b_n
    end
    return nothing
end

function cou_force_density_point!(storage::BBTMStorage, system::Peridynamics.BondSystem, 
    ::BBTMMaterial, paramhandler::Peridynamics.ParameterHandler, i::Int, mbd_m::Vector{Float64}) 

    params_i = Peridynamics.get_params(paramhandler, i)
    (; position, n_active_bonds, bond_active, temperature, b_int) = storage
    (; bonds, volume) = system

    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length 
        mof_mh =  1.0 #mbd_m[bond_id]
        Δxij = Peridynamics.get_vector_diff(position, i, j)

        params_j = Peridynamics.get_params(paramhandler, j)
        l = sqrt(Δxij[1] * Δxij[1] + Δxij[2] * Δxij[2] + Δxij[3] * Δxij[3])

        εm = (l - L) / L
        εt = 0.25 * (temperature[1, j] + temperature[1, i]) * (params_i.aph + params_j.aph)
        ε = εm - εt
        
        if ε > 0.5 * (params_i.εc + params_j.εc) && bond.fail_permit
            bond_active[bond_id] = false
        end

        n_active_bonds[i] += bond_active[bond_id]
        
        b_int[:, i] += bond_active[bond_id] * (params_i.bc + params_j.bc)/2 * ε / l * volume[j] .* Δxij * mof_mh
    end
    return nothing
end
