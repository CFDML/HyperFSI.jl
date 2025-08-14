## Here we have the modification module for thermal diffusion and coupling
## No the single mechanical part
## "new" means the "New correction ensures symmetry" for the thermal part

# New correction ensures symmetry
# and the results are basically similar
function new_modify_tm(dh::Peridynamics.AbstractDataHandler)
    ch_point_th, ch_point_m = new_cal_all_point_pd_tm(dh)
    t_factors, m_factors = new_cal_th_factor_tm(dh, ch_point_th, ch_point_m)
    return t_factors, m_factors
end

function new_cal_all_point_pd_tm(dh::Peridynamics.AbstractDataHandler)
    ch_point_th_loc = Vector{Matrix{Float64}}(undef, dh.n_chunks)
    ch_point_th = Vector{Matrix{Float64}}(undef, dh.n_chunks)
    ch_point_m_loc = Vector{Matrix{Float64}}(undef, dh.n_chunks)
    ch_point_m = Vector{Matrix{Float64}}(undef, dh.n_chunks)

    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        (; system, paramsetup) = chunk
        n_points = size(system.position, 2)
        ch_point_th_loc[chunk_id] = ones(3, chunk.system.chunk_handler.n_loc_points)
        ch_point_th[chunk_id] = ones(3, n_points)
        ch_point_m_loc[chunk_id] = ones(3, chunk.system.chunk_handler.n_loc_points)
        ch_point_m[chunk_id] = ones(3, n_points)
        tem = ones(3, n_points)
        defposition = copy(system.position)

        if paramsetup isa Peridynamics.AbstractPointParameters

            for dm in 1:3
                tem[dm, :] = 0.001.* system.position[dm, :]
                defposition[dm, :] .*= 1.001
                for i in eachindex(chunk.system.chunk_handler.loc_points)
                    pd_th = 0.0  
                    pd_me = 0.0 
                    for bond_id in system.bond_ids[i]
                        bond = system.bonds[bond_id]
                        j, L = bond.neighbor, bond.length  

                        failure = chunk.storage.bond_active[bond_id]
                        mfth = 1.0
                        mfmech = 1.0

                        for key in filter(key -> contains(string(key), "fix"), keys(chunk.psets))
                            if i in chunk.psets[key] || j in chunk.psets[key]
                                mfth = 0.0
                            end
                        end
                        
                        for key in filter(key -> contains(string(key), "tem"), keys(chunk.psets))
                            if i in chunk.psets[key] || j in chunk.psets[key]
                                mfmech = 0.0
                            end
                        end

                        pd_th +=  0.25 * failure * paramsetup.kp * (tem[dm, j] - tem[dm, i])^2 / L *
                                system.volume[j] * mfth

                        l = sqrt((defposition[1,j] - defposition[1,i])^2 + 
                                (defposition[2,j] - defposition[2,i])^2 + 
                                (defposition[3,j] - defposition[3,i])^2)

                        pd_me +=  0.25 * failure * paramsetup.bc * (l - L)^2 / L *
                                system.volume[j] * mfmech

                    end 
                    point_th = isapprox(pd_th, 0; atol=1e-6) ? 1.0 : 1e-6 * paramsetup.kc * 0.5 / pd_th              
                    ch_point_th_loc[chunk_id][dm, i] = point_th

                    if pd_me != 0
                        if abs(paramsetup.nu - 1/3) < 1e-9
                            point_m = 9/16 * 1e-6 * paramsetup.E / pd_me
                        else
                            point_m = 0.60 * 1e-6 * paramsetup.E / pd_me
                        end
                    else
                        point_m = 1.0
                    end

                    ch_point_m_loc[chunk_id][dm, i] = point_m  
                end
            end            
        else

            for dm in 1:3
                tem[dm, :] = 0.001.* system.position[dm, :]
                defposition[dm, :] .*= 1.001
                for i in eachindex(chunk.system.chunk_handler.loc_points)
                    pd_th = 0.0  
                    pd_me = 0.0 
                    params_i = Peridynamics.get_params(paramsetup, i)
                    for bond_id in system.bond_ids[i]
                        bond = system.bonds[bond_id]
                        j, L = bond.neighbor, bond.length  

                        failure = chunk.storage.bond_active[bond_id]
                        mfth = 1.0
                        mfmech = 1.0

                        params_j = Peridynamics.get_params(paramsetup, j)

                        for key in filter(key -> contains(string(key), "fix"), keys(chunk.psets))
                            if i in chunk.psets[key] || j in chunk.psets[key]
                                mfth = 0.0
                            end
                        end
                        
                        for key in filter(key -> contains(string(key), "tem"), keys(chunk.psets))
                            if i in chunk.psets[key] || j in chunk.psets[key]
                                mfmech = 0.0
                            end
                        end

                        pd_th +=  0.25 * failure * (params_i.kp + params_j.kp)/2 * (tem[dm, j] - tem[dm, i])^2 / L *
                                system.volume[j] * mfth

                        l = sqrt((defposition[1,j] - defposition[1,i])^2 + 
                                (defposition[2,j] - defposition[2,i])^2 + 
                                (defposition[3,j] - defposition[3,i])^2)

                        pd_me +=  0.25 * failure * (params_i.bc + params_j.bc)/2 * (l - L)^2 / L *
                                system.volume[j] * mfmech

                    end 
                    point_th = isapprox(pd_th, 0; atol=1e-6) ? 1.0 : 1e-6 * params_i.kc * 0.5 / pd_th              
                    ch_point_th_loc[chunk_id][dm, i] = point_th

                    if pd_me != 0
                        if abs(params_i.nu - 1/3) < 1e-6
                            point_m = 9/16 * 1e-6 * params_i.E / pd_me
                        else
                            point_m = 3/5 * 1e-6 * params_i.E / pd_me
                        end
                    else
                        point_m = 0.0
                    end

                    ch_point_m_loc[chunk_id][dm, i] = point_m  
                end
            end 
                        
        end 

    end

    all_point_th = reduce(hcat, ch_point_th_loc)
    all_point_m = reduce(hcat, ch_point_m_loc) 
    for ic in 1:dh.n_chunks
        for j in eachindex(dh.chunks[ic].system.chunk_handler.point_ids)
            ch_point_th[ic][:, j] = all_point_th[:, dh.chunks[ic].system.chunk_handler.point_ids[j]]
            ch_point_m[ic][:, j] = all_point_m[:, dh.chunks[ic].system.chunk_handler.point_ids[j]]
        end
    end
    return ch_point_th, ch_point_m
end

function new_cal_th_factor_tm(dh::Peridynamics.AbstractDataHandler, 
                            ch_point_th::Vector{Matrix{Float64}}, ch_point_m::Vector{Matrix{Float64}})
    mbd_t = Vector{Vector{Float64}}(undef, dh.n_chunks)
    mbd_m = Vector{Vector{Float64}}(undef, dh.n_chunks)
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        system = chunk.system
        point_th = ch_point_th[chunk_id]
        point_m = ch_point_m[chunk_id]
        mbd_t[chunk_id] = Vector{Float64}(undef, length(system.bonds))
        mbd_m[chunk_id] = Vector{Float64}(undef, length(system.bonds))
        for i in eachindex(chunk.system.chunk_handler.loc_points)
            for bond_id in chunk.system.bond_ids[i]
                bond = system.bonds[bond_id]
                j, L = bond.neighbor, bond.length
                Δxijx = system.position[1, j] - system.position[1, i]
                Δxijy = system.position[2, j] - system.position[2, i]
                Δxijz = system.position[3, j] - system.position[3, i]
                if abs(Δxijz) <= 1e-10
                    if abs(Δxijy) <= 1e-10
                        θ = 0.0
                    elseif abs(Δxijx) <= 1e-10
                        θ = 90 * π / 180
                    else
                        θ = atan(abs(Δxijy) / abs(Δxijx))
                    end
                    ϕ = 90 * π / 180
                    scxt = (point_th[1, i] + point_th[1, j]) / 2
                    scyt = (point_th[2, i] + point_th[2, j]) / 2
                    sczt = (point_th[3, i] + point_th[3, j]) / 2
                    scrt = sqrt(1 / ((cos(θ) * sin(ϕ))^2 / scxt^2 +
                                (sin(θ) * sin(ϕ))^2 / scyt^2 +
                                cos(ϕ)^2 / sczt^2))

                    scxm = (point_m[1, i] + point_m[1, j]) / 2
                    scym = (point_m[2, i] + point_m[2, j]) / 2
                    sczm = (point_m[3, i] + point_m[3, j]) / 2
                    scrm = sqrt(1 / ((cos(θ) * sin(ϕ))^2 / scxm^2 +
                                (sin(θ) * sin(ϕ))^2 / scym^2 +
                                cos(ϕ)^2 / sczm^2))

                elseif abs(Δxijx) <= 1e-10 && abs(Δxijy) <= 1e-10
                    sczt = (point_th[3, i] + point_th[3, j]) / 2
                    scrt = sczt

                    sczm = (point_m[3, i] + point_m[3, j]) / 2
                    scrm = sczm
                else
                    θ = atan(abs(Δxijy) / abs(Δxijx))
                    ϕ = acos(abs(Δxijz) / L)
                    scxt = (point_th[1, i] + point_th[1, j]) / 2
                    scyt = (point_th[2, i] + point_th[2, j]) / 2
                    sczt = (point_th[3, i] + point_th[3, j]) / 2
                    scrt = sqrt(1 / ((cos(θ) * sin(ϕ))^2 / scxt^2 +
                                (sin(θ) * sin(ϕ))^2 / scyt^2 +
                                cos(ϕ)^2 / sczt^2))

                    scxm = (point_m[1, i] + point_m[1, j]) / 2
                    scym = (point_m[2, i] + point_m[2, j]) / 2
                    sczm = (point_m[3, i] + point_m[3, j]) / 2
                    scrm = sqrt(1 / ((cos(θ) * sin(ϕ))^2 / scxm^2 +
                                (sin(θ) * sin(ϕ))^2 / scym^2 +
                                cos(ϕ)^2 / sczm^2))
                end
                mbd_t[chunk_id][bond_id] = scrt
                mbd_m[chunk_id][bond_id] = scrm
            end

        end        
    end
    return mbd_t, mbd_m
end


# Traditional corrections lack symmetry
# but the impact is minor
# smaller spatial dimensions can compensate, or corrections can even be omitted

function modify_t(dh::Peridynamics.AbstractDataHandler)
    ch_point_th =  cal_all_point_pd_t(dh)
    t_factors =  cal_th_factor_t(dh, ch_point_th)
    return t_factors
end

function cal_all_point_pd_t(dh::Peridynamics.AbstractDataHandler)
    ch_point_th_loc = Vector{Vector{Float64}}(undef, dh.n_chunks)
    ch_point_th = Vector{Vector{Float64}}(undef, dh.n_chunks)

    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        (; system, paramsetup) = chunk
        n_points = size(system.position, 2)
        ch_point_th_loc[chunk_id] =Vector{Float64}(undef, chunk.system.chunk_handler.n_loc_points)
        ch_point_th[chunk_id] = Vector{Float64}(undef, n_points)

        tem =  0.001.*[sum(system.position[:, i]) for i in 1:n_points]

        if paramsetup isa Peridynamics.AbstractPointParameters
            for i in eachindex(chunk.system.chunk_handler.loc_points)
                pd_th = 0.0   
                for bond_id in system.bond_ids[i]
                    bond = system.bonds[bond_id]
                    j, L = bond.neighbor, bond.length  
                    pd_th +=  0.25 * paramsetup.kp * (tem[j] - tem[i])^2 / L *
                            system.volume[j] 
                end 
                point_th = isapprox(paramsetup.nu, 1/3; atol=1e-5) ? 1e-6 * paramsetup.kc * 1.0 / pd_th  : 1e-6 * paramsetup.kc * 1.5 / pd_th              
                ch_point_th_loc[chunk_id][i] = point_th
            end            
        else
            for i in eachindex(chunk.system.chunk_handler.loc_points)
                pd_th = 0.0   
                params_i = Peridynamics.get_params(paramsetup, i)
                for bond_id in system.bond_ids[i]
                    bond = system.bonds[bond_id]
                    j, L = bond.neighbor, bond.length  
                    params_j = Peridynamics.get_params(paramsetup, j)
                    pd_th +=  0.25 * (params_i.kp + params_j.kp)/2 * (tem[j] - tem[i])^2 / L *
                            system.volume[j] 
                end 
                point_th = isapprox(params_i.nu, 1/3; atol=1e-5) ? 1e-6 * params_i.kc * 1.0 / pd_th  : 1e-6 * params_i.kc * 1.5 / pd_th              
                ch_point_th_loc[chunk_id][i] = point_th
            end
        end

    end

    all_point_th = vcat(ch_point_th_loc...)
    for ic in 1:dh.n_chunks
        for j in eachindex(dh.chunks[ic].system.chunk_handler.point_ids)
            ch_point_th[ic][j] = all_point_th[dh.chunks[ic].system.chunk_handler.point_ids[j]]
        end
    end
    return ch_point_th
end

function cal_th_factor_t(dh::Peridynamics.AbstractDataHandler, ch_point_th::Vector{Vector{Float64}})
    mbd_t = Vector{Vector{Float64}}(undef, dh.n_chunks)
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        system = chunk.system
        point_th = ch_point_th[chunk_id]
        mbd_t[chunk_id] = Vector{Float64}(undef, length(system.bonds))
        for i in eachindex(chunk.system.chunk_handler.loc_points)
            for bond_id in chunk.system.bond_ids[i]
                bond = chunk.system.bonds[bond_id]
                j, L = bond.neighbor, bond.length
                mbd_t[chunk_id][bond_id] = 0.5 * (point_th[i] + point_th[j])
            end
        end
    end
    return mbd_t
end

# New correction ensures symmetry
# and the results are basically similar
function new_modify_t(dh::Peridynamics.AbstractDataHandler)
    ch_point_th = new_cal_all_point_pd_t(dh)
    t_factors = new_cal_th_factor_t(dh, ch_point_th)
    return t_factors
end

function new_cal_all_point_pd_t(dh::Peridynamics.AbstractDataHandler)
    ch_point_th_loc = Vector{Matrix{Float64}}(undef, dh.n_chunks)
    ch_point_th = Vector{Matrix{Float64}}(undef, dh.n_chunks)

    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        (; system, paramsetup) = chunk
        n_points = size(system.position, 2)
        ch_point_th_loc[chunk_id] = ones(3, chunk.system.chunk_handler.n_loc_points)
        ch_point_th[chunk_id] = ones(3, n_points)
        tem = ones(3, n_points)

        for dm in 1:3
            tem[dm, :] = 0.001.* system.position[dm, :]

            if paramsetup isa Peridynamics.AbstractPointParameters
                for i in eachindex(chunk.system.chunk_handler.loc_points)
                    pd_th = 0.0   
                    for bond_id in system.bond_ids[i]
                        bond = system.bonds[bond_id]
                        j, L = bond.neighbor, bond.length  
                        pd_th +=  0.25 * paramsetup.kp * (tem[dm, j] - tem[dm, i])^2 / L *
                                system.volume[j] 
                    end 
                    point_th = isapprox(pd_th, 0; atol=1e-8) ? 1.0 : 1e-6 * paramsetup.kc * 0.5 / pd_th             
                    ch_point_th_loc[chunk_id][dm, i] = point_th
                end            
            else
                for i in eachindex(chunk.system.chunk_handler.loc_points)
                    pd_th = 0.0   
                    params_i = Peridynamics.get_params(paramsetup, i)
                    for bond_id in system.bond_ids[i]
                        bond = system.bonds[bond_id]
                        j, L = bond.neighbor, bond.length  
                        params_j = Peridynamics.get_params(paramsetup, j)
                        pd_th +=  0.25 * (params_i.kp + params_j.kp)/2 * (tem[dm, j] - tem[dm, i])^2 / L *
                                system.volume[j] 
                    end 
                    point_th = isapprox(pd_th, 0; atol=1e-8) ? 1.0 : 1e-6 * params_i.kc * 0.5 / pd_th           
                    ch_point_th_loc[chunk_id][dm, i] = point_th
                end
            end            
        end
    end

    all_point_th = reduce(hcat, ch_point_th_loc)
    for ic in 1:dh.n_chunks
        for j in eachindex(dh.chunks[ic].system.chunk_handler.point_ids)
            ch_point_th[ic][:, j] = all_point_th[:, dh.chunks[ic].system.chunk_handler.point_ids[j]]
        end
    end
    return ch_point_th
end

function new_cal_th_factor_t(dh::Peridynamics.AbstractDataHandler, ch_point_th::Vector{Matrix{Float64}})
    mbd_t = Vector{Vector{Float64}}(undef, dh.n_chunks)
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        system = chunk.system
        point_th = ch_point_th[chunk_id]
        mbd_t[chunk_id] = Vector{Float64}(undef, length(system.bonds))
        for i in eachindex(chunk.system.chunk_handler.loc_points)
            for bond_id in chunk.system.bond_ids[i]
                bond = system.bonds[bond_id]
                j, L = bond.neighbor, bond.length
                Δxijx = system.position[1, j] - system.position[1, i]
                Δxijy = system.position[2, j] - system.position[2, i]
                Δxijz = system.position[3, j] - system.position[3, i]
                if abs(Δxijz) <= 1e-10
                    if abs(Δxijy) <= 1e-10
                        θ = 0.0
                    elseif abs(Δxijx) <= 1e-10
                        θ = 90 * π / 180
                    else
                        θ = atan(abs(Δxijy) / abs(Δxijx))
                    end
                    ϕ = 90 * π / 180
                    scx = (point_th[1, i] + point_th[1, j]) / 2
                    scy = (point_th[2, i] + point_th[2, j]) / 2
                    scz = (point_th[3, i] + point_th[3, j]) / 2
                    scr = sqrt(1 / ((cos(θ) * sin(ϕ))^2 / scx^2 +
                                (sin(θ) * sin(ϕ))^2 / scy^2 +
                                cos(ϕ)^2 / scz^2))
                elseif abs(Δxijx) <= 1e-10 && abs(Δxijy) <= 1e-10
                    scz = (point_th[3, i] + point_th[3, j]) / 2
                    scr = scz
                else
                    θ = atan(abs(Δxijy) / abs(Δxijx))
                    ϕ = acos(abs(Δxijz) / L)
                    scx = (point_th[1, i] + point_th[1, j]) / 2
                    scy = (point_th[2, i] + point_th[2, j]) / 2
                    scz = (point_th[3, i] + point_th[3, j]) / 2
                    scr = sqrt(1 / ((cos(θ) * sin(ϕ))^2 / scx^2 +
                                (sin(θ) * sin(ϕ))^2 / scy^2 +
                                cos(ϕ)^2 / scz^2))
                end
                mbd_t[chunk_id][bond_id] = scr
            end
        end        
    end
    return mbd_t
end


