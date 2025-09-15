"""
Bond-Based thermomechanics Module with Temperature-Dependent Properties and Ablation

This module extends the bond-based thermomechanical model 
by incorporating Ablation Behavior.

··This module provides functionality for ablation, primarily for material point removal.
··Uses an 'exist' state to indicate ablation status, which relates to ablation amount.
··Specific chemical processes are configured in the ablation file, with flexible complexity.
··Core functionality identifies new boundary points, 
  and more critically, establishes new watertight boundaries for subsequent fluid-structure interaction.
··Ablated points can be modeled to describe subsequent effects.
··The variation of properties with temperature is under consideration.

··BBTMA means Bond-Based thermomechanics with Temperature-Dependent Properties and Ablation.
"""
@inline thermomech_ablation_kwargs() = (:E, :nu, :G, :K, :lambda, :mu, :rho, :horizon, 
:Gc, :epsilon_c, :kc, :aph, :cv, :rft, :h, :hσ, :hϵ, :tem∞, :thick, 
:kp_T, :cv_T, :aph_T, :ρ_T, :bc_T, :ablation_v, :energy)

struct BBTMAMaterial{Correction,DM} <: Peridynamics.AbstractBondSystemMaterial{Correction}
    dmgmodel::DM
    function BBTMAMaterial{C}(dmgmodel::DM) where {C,DM}
        new{C,DM}(dmgmodel)
    end
end

function BBTMAMaterial{C}(; dmgmodel::Peridynamics.AbstractDamageModel=CriticalStretch()) where {C}
    return BBTMAMaterial{C}(dmgmodel)
end
BBTMAMaterial(; kwargs...) = BBTMAMaterial{NoCorrection}(; kwargs...)

struct BBTMAPointParameters <: Peridynamics.AbstractPointParameters
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
    kp_T::Function # function describing PD thermal_conductity variation with temperature f(kp,T)
    cv_T::Function # function describing cv variation with temperature f(cv,T)
    aph_T::Function # function describing aph variation with temperature f(aph,T)
    ρ_T::Function # function describing density variation with temperature f(ρ,T)
    bc_T::Function # function describing bond constant variation with temperature f(bc,T)
    ablation_v::Function # function describing ablation rate f(vol_reference, T)
    energy::Float64 # total energy release from ablation
end

function const_ρ(rho, T)
    return rho
end

function const_bc(E, bc, T)
    return bc
end

function BBTMAPointParameters(mat::BBTMAMaterial, p::Dict{Symbol,Any})
    par = Peridynamics.get_given_elastic_params(p)
    (; E, nu, G, K, λ, μ) = par

    if haskey(p, :thick)
        p[:nu] = 1/3
    else 
        p[:nu] = 1/4    
    end
   
    (; δ, rho, E, nu, G, K, λ, μ) = Peridynamics.get_required_point_parameters(mat, p)
    Gc, εc = Peridynamics.get_frac_params(mat.dmgmodel, p, δ, K)
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
    ρ_T = get(p, :ρ_T, const_ρ)
    bc_T = get(p, :bc_T, const_bc)
    ablation_v = get(p, :ablation_v, const_ablation_v)

    return BBTMAPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc, kc, kp, aph, cv, rft, h, hσ, hϵ, tem∞,
        kp_T, cv_T, aph_T, ρ_T, bc_T, ablation_v, energy)
end

@Peridynamics.params BBTMAMaterial BBTMAPointParameters

@Peridynamics.storage BBTMAMaterial struct BBTMAStorage <: Peridynamics.AbstractStorage
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
    @pointfield ablation_deep::Matrix{Float64}
    @lthfield exist::Matrix{Int} # 1: exist, 0: ablated
end


function Peridynamics.init_field_solver(::Thermomechstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:position})
    return copy(system.position)
end

function Peridynamics.init_field_solver(::Thermomechstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:displacement})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:velocity})
    return zeros(3, size(system.position, 2))
end

function Peridynamics.init_field_solver(::Thermomechstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:velocity_half})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:acceleration})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:b_int})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:b_ext})
    return zeros(3, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:temperature})
    return zeros(1, size(system.position, 2))
end

function Peridynamics.init_field_solver(::Thermomechstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:pflux})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:hsource})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:ablation_deep})
    return zeros(1, Peridynamics.get_n_loc_points(system))
end

function Peridynamics.init_field_solver(::Thermomechstep_ablation, system::Peridynamics.AbstractSystem, ::Val{:exist})
    return ones(Int, 1, size(system.position, 2))
end

function Peridynamics.Peridynamics.init_field(::BBTMAMaterial, ::Peridynamics.AbstractTimeSolver, system::Peridynamics.BondSystem,
                    ::Val{:bond_stretch})
    return zeros(Peridynamics.get_n_bonds(system))
end

function Peridynamics.Peridynamics.allowed_material_kwargs(::BBTMAMaterial)
    return (thermomech_ablation_kwargs())
end

function new_calc_force_density_ablation!(dh::Peridynamics.ThreadsBodyDataHandler, t, Δt)
    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_loc_to_halo!(dh, chunk_id)
        new_calc_force_density_ablation!(dh.chunks[chunk_id], t, Δt)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_halo_to_loc!(dh, chunk_id)
    end
    return nothing
end

function new_calc_force_density_ablation!(chunk::Peridynamics.AbstractBodyChunk{<:Peridynamics.AbstractBondSystem}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    (; dmgmodel) = mat
    storage.b_int .= 0
    storage.n_active_bonds .= 0
    for point_id in Peridynamics.each_point_idx(chunk)
        cou_calc_failure!(storage, system, mat, dmgmodel, paramsetup, point_id)
        Peridynamics.calc_damage!(storage, system, mat, dmgmodel, paramsetup, point_id)
        new_force_density_point!(storage, system, mat, paramsetup, point_id)
    end
    Peridynamics.nancheck(chunk, t, Δt)
    return nothing
end

function cou_calc_failure!(storage::BBTMAStorage, system::Peridynamics.BondSystem,
                       ::BBTMAMaterial, ::Peridynamics.CriticalStretch,
                       paramsetup::BBTMAPointParameters, i)
    (; εc, aph) = Peridynamics.get_params(paramsetup, i)
    (; position, n_active_bonds, bond_active, bond_stretch) = storage
    (; bonds) = system
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        if storage.exist[1, i] == 0 || storage.exist[1, j] == 0
            bond_active[bond_id] = false
            bond_stretch[bond_id] = 0.0
        else
            Δxij = Peridynamics.get_vector_diff(position, i, j)
            l = norm(Δxij)
            ε = (l - L) / L - 0.5 * (storage.temperature[1, j] + storage.temperature[1, i]) * aph
            bond_stretch[bond_id] = ε / l # note that this is  ε / l!
            if ε > εc && bond.fail_permit
                bond_active[bond_id] = false
            end
        end
        n_active_bonds[i] += bond_active[bond_id]
    end
    return nothing
end

function cou_calc_failure!(storage::BBTMAStorage, system::Peridynamics.BondSystem,
                       ::BBTMAMaterial, ::Peridynamics.CriticalStretch,
                       paramhandler::Peridynamics.ParameterHandler, i)
    params_i = Peridynamics.get_params(paramhandler, i)
    (; position, n_active_bonds, bond_active, bond_stretch) = storage
    (; bonds) = system
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        if storage.exist[1, i] == 0 || storage.exist[1, j] == 0
            bond_active[bond_id] = false
            bond_stretch[bond_id] = 0.0
        else
            params_j = Peridynamics.get_params(paramhandler, j)
            Δxij = Peridynamics.get_vector_diff(position, i, j)
            l = norm(Δxij)
            ε = (l - L) / L - 0.5 * (storage.temperature[1, j] + storage.temperature[1, i]) * (params_i.aph + params_j.aph) / 2
            bond_stretch[bond_id] = ε / l # note that this is  ε / l!
            εcm = min(params_i.εc, params_j.εc)^2 / max(params_i.εc, params_j.εc) # interfacial critical stretch
            if ε > εcm && bond.fail_permit
                bond_active[bond_id] = false
            end
        end
        n_active_bonds[i] += bond_active[bond_id]
    end
    return nothing
end

function new_force_density_point!(storage::BBTMAStorage, system::Peridynamics.BondSystem, ::BBTMAMaterial,
                              params::BBTMAPointParameters, i)
    (; position, bond_stretch, bond_active, b_int) = storage
    (; bonds, volume) = system
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        Δxij = Peridynamics.get_vector_diff(position, i, j)
        ε = bond_stretch[bond_id]
        ω = bond_active[bond_id] 
        b = ω * params.bc * ε * volume[j] .* Δxij
        Peridynamics.update_add_vector!(b_int, i, b)
    end
    return nothing
end

function new_force_density_point!(storage::BBTMAStorage, system::Peridynamics.BondSystem, ::BBTMAMaterial,
                              paramhandler::Peridynamics.ParameterHandler, i)
    (; position, bond_stretch, bond_active, b_int) = storage
    (; bonds, volume) = system
    params_i = Peridynamics.get_params(paramhandler, i)
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        Δxij = Peridynamics.get_vector_diff(position, i, j)
        ε = bond_stretch[bond_id]
        params_j = Peridynamics.get_params(paramhandler, j)
        ω = bond_active[bond_id] 
        b = ω * (params_i.bc + params_j.bc) / 2 * ε * volume[j] .* Δxij
        Peridynamics.update_add_vector!(b_int, i, b)
    end
    return nothing
end

function cou_calc_pflux!(chunk::Peridynamics.AbstractBodyChunk)
    (; system, mat, paramsetup, storage) = chunk
    storage.pflux .= 0.0
    for point_id in eachindex(chunk.system.chunk_handler.loc_points)
        cou_pflux_point!(storage, system, mat, paramsetup, point_id)
    end
    return nothing
end

function cou_pflux_point!(storage::BBTMAStorage, system::Peridynamics.BondSystem, 
    ::BBTMAMaterial, param::BBTMAPointParameters, i::Int) 

    if storage.exist[1, i] == 0
        storage.pflux[1, i] = 0.0
    else
        for bond_id in system.bond_ids[i]
            bond = system.bonds[bond_id]
            j, L = bond.neighbor, bond.length

            if storage.exist[1, j] == 0
                Δflux = 0.0
            else
                Δxijx = storage.position[1, j] - storage.position[1, i]
                Δxijy = storage.position[2, j] - storage.position[2, i]
                Δxijz = storage.position[3, j] - storage.position[3, i]
                Δvijx = storage.velocity[1, j] - storage.velocity[1, i]
                Δvijy = storage.velocity[2, j] - storage.velocity[2, i]
                Δvijz = storage.velocity[3, j] - storage.velocity[3, i]

                l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)
                ev = (Δxijx * Δvijx + Δxijy * Δvijy + Δxijz * Δvijz) / l
                Δtem = storage.temperature[1, j] - storage.temperature[1, i]

                kp_tm = 0.5*(param.kp_T(param.kc, param.kp, storage.temperature[1, i]) 
                        + param.kp_T(param.kc, param.kp, storage.temperature[1, j]))

                bc_tm = 0.5*(param.bc_T(param.E, param.bc, storage.temperature[1, i])
                        + param.bc_T(param.E, param.bc, storage.temperature[1, j]))

                aph_tm = 0.5*(param.aph_T(param.aph, storage.temperature[1, i])
                        + param.aph_T(param.aph, storage.temperature[1, j]))

                Δflux = storage.bond_active[bond_id] * (kp_tm * Δtem / L -
                        ev * 0.5 * param.rft * bc_tm * aph_tm) * system.volume[j]
            end

            storage.pflux[1, i] += Δflux
        end
    end
    return nothing
end

function cou_pflux_point!(storage::BBTMAStorage, system::Peridynamics.BondSystem, 
    ::BBTMAMaterial, paramhandler::Peridynamics.ParameterHandler, i::Int) 

    if storage.exist[1, i] == 0
        storage.pflux[1, i] = 0.0
    else
        params_i = Peridynamics.get_params(paramhandler, i)
        for bond_id in system.bond_ids[i]
            bond = system.bonds[bond_id]
            j, L = bond.neighbor, bond.length
            
            if storage.exist[1, j] == 0
                Δflux = 0.0
            else
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

                kp_tm = 0.5*(params_i.kp_T(params_i.kc, params_i.kp, storage.temperature[1, i]) 
                        + params_j.kp_T(params_j.kc, params_j.kp, storage.temperature[1, j]))


                bc_tm = 0.5*(params_i.bc_T(params_i.E, params_i.bc, storage.temperature[1, i])
                        + params_j.bc_T(params_j.E, params_j.bc, storage.temperature[1, j]))

                aph_tm = 0.5*(params_i.aph_T(params_i.aph, storage.temperature[1, i])
                        + params_j.aph_T(params_j.aph, storage.temperature[1, j]))

                Δflux = storage.bond_active[bond_id] * (kp_tm * Δtem / L - ev * 0.5 * 0.5 * (params_i.rft + params_j.rft) 
                                * bc_tm *  aph_tm) * system.volume[j] 
            end
            storage.pflux[1, i] += Δflux
        end
    end
    return nothing
end

function cou_force_density_point!(storage::BBTMAStorage, system::Peridynamics.BondSystem, 
    ::BBTMAMaterial, param::BBTMAPointParameters, i::Int) 

    if storage.exist[1, i] == 0
        storage.b_int[:, i] .= 0.0
    else

        (; position, bond_stretch, bond_active, b_int) = storage
        (; bonds, volume) = system

        for bond_id in Peridynamics.each_bond_idx(system, i)
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length 
            bc_tm = 0.5*(param.bc_T(param.E, param.bc, storage.temperature[1, i])
                        + param.bc_T(param.E, param.bc, storage.temperature[1, j]))

            if storage.exist[1, j] == 0
                b_n = [0.0; 0.0; 0.0]
            else
                Δxij = Peridynamics.get_vector_diff(position, i, j)          
                b_n = bond_active[bond_id] * bc_tm * bond_stretch[bond_id] * volume[j] .* Δxij
            end
            b_int[:, i] += b_n
        end
    end
    return nothing
end

function cou_force_density_point!(storage::BBTMAStorage, system::Peridynamics.BondSystem, 
    ::BBTMAMaterial, paramhandler::Peridynamics.ParameterHandler, i::Int) 

    if storage.exist[1, i] == 0
        storage.b_int[:, i] .= 0.0
    else
        params_i = Peridynamics.get_params(paramhandler, i)
        (; position, bond_stretch, bond_active, b_int) = storage
        (; bonds, volume) = system

        for bond_id in Peridynamics.each_bond_idx(system, i)
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length 

            if storage.exist[1, j] == 0
                b_n = [0.0; 0.0; 0.0]
            else
                Δxij = Peridynamics.get_vector_diff(position, i, j)
                params_j = Peridynamics.get_params(paramhandler, j)

                bc_tm = 0.5*(params_i.bc_T(params_i.E, params_i.bc, storage.temperature[1, i])
                        + params_j.bc_T(params_j.E, params_j.bc, storage.temperature[1, j]))

                b_n = bond_active[bond_id] * bc_tm * bond_stretch[bond_id] * volume[j] .* Δxij
            end

            b_int[:, i] += b_n
        end
    end
    return nothing
end
