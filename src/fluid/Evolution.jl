function ibm_evolve!(
    KS::SolverSet,
    ctr::KB.AM{TC},
    a1face::KB.AM{TF},
    a2face::KB.AM{TF},
    flags,
    ghost_ids,
    ghost_ctr,
    dt;
    mode = KB.symbolize(KS.set.flux)::Symbol,
    bc = KB.symbolize(KS.set.boundary),
) where {TC<:Union{ControlVolume,ControlVolume2D},TF<:Union{Interface,Interface2D}}

    if mode == :gks
        ibm_gks_evolve!(KS, ctr, a1face, a2face,flags, ghost_ids, ghost_ctr, dt)
    elseif mode == :hllc
        ibm_hllc_evolve!(KS, ctr, a1face, a2face,flags, ghost_ids, ghost_ctr, dt)
    else
        error("其他求解器正在开发中")
    end

    return nothing
end

function ibm_gks_evolve!(
    KS::SolverSet,
    ctr::KB.AM{TC},
    a1face::KB.AM{TF},
    a2face::KB.AM{TF},
    flags,
    ghost_ids,
    ghost_ctr,
    dt
) where {TC<:Union{ControlVolume,ControlVolume2D},TF<:Union{Interface,Interface2D}}
    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    if firstindex(KS.ps.x[:, 1]) < 1
        idx0 = 1
        idx1 = nx + 1
    else
        idx0 = 2
        idx1 = nx
    end
    if firstindex(KS.ps.y[1, :]) < 1
        idy0 = 1
        idy1 = ny + 1
    else
        idy0 = 2
        idy1 = ny
    end    

    # x direction
    @inbounds @threads for j = 1:ny
        for i = idx0:idx1
            dxL = 0.5 * dx[i-1, j]
            dxR = 0.5 * dx[i, j]
            n = KS.ps.n[i-1, j, 2]
            len = KS.ps.areas[i-1, j, 2]
            axis = 1
            swL = ctr[i-1, j].sw[:, axis]
            swR = ctr[i, j].sw[:, axis]

            if flags[i-1, j] == -2 && flags[i, j] == 1
                direc = CartesianIndex(1, 0)
                indices = findall(p -> p == CartesianIndex(i-1, j), ghost_ids)
                @assert size(indices, 1) == 1
                indice = indices[1]
                wL = ghost_ctr[indice][direc].w .+ dxL .* swL
                wR = ctr[i, j].w .- dxR .* swR
            elseif flags[i-1, j] == 1 && flags[i, j] == -2
                direc = CartesianIndex(-1, 0)
                indices = findall(p -> p == CartesianIndex(i, j), ghost_ids)
                @assert size(indices, 1) == 1
                indice = indices[1]
                wL = ctr[i-1, j].w
                wR = ghost_ctr[indice][direc].w
            else
                wL = ctr[i-1, j].w
                wR = ctr[i, j].w
            end

            flux_gks!(
                a1face[i, j].fw,
                local_frame(wL .+ dxL .* swL, n),
                local_frame(wR .- dxR .* swR, n),
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                dt,
                dxL,
                dxR,
                len,
                swL,
                swR
            )

            a1face[i, j].fw .= global_frame(a1face[i, j].fw, n)
        end
    end

    # y direction
    @inbounds @threads for j = idy0:idy1
        for i = 1:nx
            dxL = 0.5 * dy[i, j-1]
            dxR = 0.5 * dy[i, j]
            n = KS.ps.n[i, j-1, 3]
            len = KS.ps.areas[i, j-1, 3]
            axis = 2
            swL = ctr[i, j-1].sw[:, axis]
            swR = ctr[i, j].sw[:, axis]

            if flags[i, j-1] == -2 && flags[i, j] == 1
                direc = CartesianIndex(0, 1)
                indices = findall(p -> p == CartesianIndex(i, j-1), ghost_ids)
                @assert size(indices, 1) == 1
                indice = indices[1]
                wL = ghost_ctr[indice][direc].w
                wR = ctr[i, j].w
            elseif flags[i, j-1] == 1 && flags[i, j] == -2
                direc = CartesianIndex(0, -1)
                indices = findall(p -> p == CartesianIndex(i, j), ghost_ids)
                @assert size(indices, 1) == 1
                indice = indices[1]
                wL = ctr[i, j-1].w
                wR = ghost_ctr[indice][direc].w
            else
                wL = ctr[i, j-1].w
                wR = ctr[i, j].w
            end

            flux_gks!(
                a2face[i, j].fw,
                local_frame(wL .+ dxL .* swL, n),
                local_frame(wR .- dxR .* swR, n),
                KS.gas.K,
                KS.gas.γ,
                KS.gas.μᵣ,
                KS.gas.ω,
                dt,
                dxL,
                dxR,
                len,
                swL,
                swR
            )

            a2face[i, j].fw .= global_frame(a2face[i, j].fw, n)
        end
    end
end

function ibm_hllc_evolve!(
    KS::SolverSet,
    ctr::KB.AM{TC},
    a1face::KB.AM{TF},
    a2face::KB.AM{TF},
    flags,
    ghost_ids,
    ghost_ctr,
    dt
) where {TC<:Union{ControlVolume,ControlVolume2D},TF<:Union{Interface,Interface2D}}
    nx, ny, dx, dy = begin
        if KS.ps isa CSpace2D
            KS.ps.nr, KS.ps.nθ, KS.ps.dr, KS.ps.darc
        else
            KS.ps.nx, KS.ps.ny, KS.ps.dx, KS.ps.dy
        end
    end

    if firstindex(KS.ps.x[:, 1]) < 1
        idx0 = 1
        idx1 = nx + 1
    else
        idx0 = 2
        idx1 = nx
    end
    if firstindex(KS.ps.y[1, :]) < 1
        idy0 = 1
        idy1 = ny + 1
    else
        idy0 = 2
        idy1 = ny
    end

    # x direction
    @inbounds @threads for j = 1:ny
        for i = idx0:idx1
            n = KS.ps.n[i-1, j, 2]
            len = KS.ps.areas[i-1, j, 2]

            if flags[i-1, j] == -2 && flags[i, j] == 1
                direc = CartesianIndex(1, 0)
                indices = findall(p -> p == CartesianIndex(i-1, j), ghost_ids)
                @assert size(indices, 1) == 1
                indice = indices[1]
                KitBase.flux_hllc!(
                    KS,
                    a1face[i, j],
                    ghost_ctr[indice][direc],
                    ctr[i, j],
                    (0.5 .* dx[i-1, j], 0.5 .* dx[i, j], len, n, 1),
                    dt,
                )
            elseif flags[i-1, j] == 1 && flags[i, j] == -2
                direc = CartesianIndex(-1, 0)
                indices = findall(p -> p == CartesianIndex(i, j), ghost_ids)
                @assert size(indices, 1) == 1
                indice = indices[1]
                KitBase.flux_hllc!(
                    KS,
                    a1face[i, j],
                    ctr[i-1, j],
                    ghost_ctr[indice][direc],
                    (0.5 .* dx[i-1, j], 0.5 .* dx[i, j], len, n, 1),
                    dt,
                )
            else
                KitBase.flux_hllc!(
                    KS,
                    a1face[i, j],
                    ctr[i-1, j],
                    ctr[i, j],
                    (0.5 .* dx[i-1, j], 0.5 .* dx[i, j], len, n, 1),
                    dt,
                )
            end
        end
    end

    # y direction
    @inbounds @threads for j = idy0:idy1
        for i = 1:nx
            n = KS.ps.n[i, j-1, 3]
            len = KS.ps.areas[i, j-1, 3]

            if flags[i, j-1] == -2 && flags[i, j] == 1
                direc = CartesianIndex(0, 1)
                indices = findall(p -> p == CartesianIndex(i, j-1), ghost_ids)
                @assert size(indices, 1) == 1
                indice = indices[1]
                KitBase.flux_hllc!(
                    KS,
                    a2face[i, j],
                    ghost_ctr[indice][direc],
                    ctr[i, j],
                    (0.5 .* dy[i, j-1], 0.5 .* dy[i, j], len, n, 2),
                    dt,
                )
            elseif flags[i, j-1] == 1 && flags[i, j] == -2
                direc = CartesianIndex(0, -1)
                indices = findall(p -> p == CartesianIndex(i, j), ghost_ids)
                @assert size(indices, 1) == 1
                indice = indices[1]
                KitBase.flux_hllc!(
                    KS,
                    a2face[i, j],
                    ctr[i, j-1],
                    ghost_ctr[indice][direc],
                    (0.5 .* dy[i, j-1], 0.5 .* dy[i, j], len, n, 2),
                    dt,
                )
            else
                KitBase.flux_hllc!(
                    KS,
                    a2face[i, j],
                    ctr[i, j-1],
                    ctr[i, j],
                    (0.5 .* dy[i, j-1], 0.5 .* dy[i, j], len, n, 2),
                    dt,
                )
            end
        end
    end
end