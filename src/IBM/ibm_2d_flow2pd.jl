function update_pressure_on_edge!(ps::KB.AbstractPhysicalSpace2D, ctr::KB.AM, 
    pos::Matrix{Float64}, area::Vector{Float64}, press_3bc::Matrix{Float64}, ib::IBM2D)
    flags = ib.flags
    edges, nodes, idpds, xbis, nbis, xips, ip_nids, ip_bids = 
        ib.f2pd.edges, ib.f2pd.nodes, ib.f2pd.idpd, ib.f2pd.xbi, ib.f2pd.nbi, ib.f2pd.xip, ib.f2pd.idn, ib.f2pd.idb

    p_ref = ib.rv.pressure
    
    for i in axes(press_3bc, 2)
        press_3bc[:, i] .= [0.0, 0.0, 0.0]
    end

    @threads for iter in axes(edges, 1)
        edge = edges[iter]
        x, y = xips[iter]
        if !(ps.x0 < x < ps.x1 && ps.y0 < y < ps.y1)
            continue
        end

        node = nodes[iter]
        pf = [1, x, y, x * y]
        direc = -1.0 .* xbis[iter]
        nids = ip_nids[iter]
        bids = ip_bids[iter]

        if size(nids, 1) > 1 && size(nids, 1) + size(bids, 1) == 4
            w1 = [0.5 * ctr[idx].prim[1] / ctr[idx].prim[4] for idx in nids]
            w = w1

            C = bilinear_coeffs(ps, edge, nids, bids, w)
            P1 = C' * pf
        elseif size(nids, 1) >= 1
            idx1 = nids[1]
            ρ1, U1, V1, λ1 = ctr[idx1].prim
            P1 = 0.5 * ρ1 / λ1
        else
            error("No fluid nodes found around the IP node!")
        end

        press_3bc[1:2, node] .+= (p_ref * P1 * norm(edge[1] .- edge[2]) / area[node]) * direc 
    end
end

function update_hsource_on_edge!(ps::KB.AbstractPhysicalSpace2D, ctr::KB.AM, 
    pos::Matrix{Float64}, area::Vector{Float64}, temperature::Vector{Float64}, hsource_1bc::Matrix{Float64}, ib::IBM2D)
    flags = ib.flags
    edges, nodes, pd_ids, xbis, nbis, xips, ip_nids, ip_bids = 
        ib.f2pd.edges, ib.f2pd.nodes, ib.f2pd.idpd, ib.f2pd.xbi, ib.f2pd.nbi, ib.f2pd.xip, ib.f2pd.idn, ib.f2pd.idb

    tem_ref = ib.rv.temperature
    l_ref = ib.rv.length
    k = ib.rv.k
    
    for i in axes(hsource_1bc, 2)
        hsource_1bc[i] = 0.0
    end

    @threads for iter in axes(edges, 1)
        edge = edges[iter]
        x, y = xips[iter]
        node = nodes[iter]
        pd_id = pd_ids[iter]
        pf = [1, x, y, x * y]
        nx, ny = -1.0 .* xbis[iter]
        nids = ip_nids[iter]
        bids = ip_bids[iter]
        t1 = 0.5 * (temperature[pd_id[1][1]] + temperature[pd_id[1][2]]) / tem_ref 
        t2 = 0.5 * (temperature[pd_id[2][1]] + temperature[pd_id[2][2]]) / tem_ref

        if size(nids, 1) > 1 && size(nids, 1) + size(bids, 1) == 4
            w1 = [1.0 / ctr[idx].prim[4] for idx in nids]
            w2 = ones(size(bids, 1))
            for i in axes(bids, 1)
                xb, yb = foot([ps.x[bids[i]], ps.y[bids[i]]], edge[1], edge[2])
                w2[i] = temperature_interpolation([xb, yb], edge[1], t1, edge[2], t2)
            end
            w = [w1; w2]

            C = bilinear_coeffs(ps, edge, nids, bids, w)
            hs = (C[2] * nx + C[3] * ny + C[4] * nx * y + C[4] * ny * x)
        elseif size(nids, 1) >= 1
            idx1 = nids[1]
            ρ1, U1, V1, λ1 = ctr[idx1].prim
            T1 = 1 / λ1
            xb, yb = foot([ps.x[idx1], ps.y[idx1]], edge[1], edge[2])
            T2 = temperature_interpolation([xb, yb], edge[1], t1, edge[2], t2)
            hs = (T2 - T1) / norm([xb, yb] .- [ps.x[idx1], ps.y[idx1]])
        else
            error("No fluid nodes found around the IP node!")
        end

        hsource_1bc[node] -= (hs * norm(edge[1] .- edge[2]) / area[node]) * (k * tem_ref / l_ref)
    end
end

function bilinear_coeffs(ps::KB.AbstractPhysicalSpace2D, xbms::Vector{Vector{Float64}}, nids::Vector{CartesianIndex}, bids::Vector{CartesianIndex}, W0::Vector{Float64})
    M = Matrix{Float64}[]
    for id in nids
        xf = ps.x[id]
        yf = ps.y[id]
        push!(M, [1 xf yf xf * yf])
    end

    if size(W0, 1) > size(nids, 1)
        for id in bids
            xg = ps.x[id]
            yg = ps.y[id]
            xb, yb = foot([xg, yg], xbms[1], xbms[2])
            push!(M, [1 xb yb xb * yb])
        end
        W = W0
    else
        for id in bids
            xg = ps.x[id]
            yg = ps.y[id]
            xb, yb = foot([xg, yg], xbms[1], xbms[2])
            nx, ny = ([xb, yb] .- [xg, yg]) / norm([xb, yb] .- [xg, yg])
            push!(M, [0 nx ny xb * ny + yb * nx])
        end
        W = [W0; zeros(size(bids, 1))]
    end
    M = vcat(M...)

    return M \ W
end

function temperature_interpolation(p::Vector{Float64}, p1::Vector{Float64}, t1::Float64, p2::Vector{Float64}, t2::Float64)
    x, y = p
    x1, y1 = p1
    x2, y2 = p2
    
    # Calculate squared length of segment p1p2
    L2 = (x2 - x1)^2 + (y2 - y1)^2
    
    # Calculate vectors from p to p1 and p2
    v1 = [x - x1, y - y1]
    v2 = [x2 - x1, y2 - y1]
    
    # Calculate parameter t
    t = dot(v1, v2) / L2
    
    # Interpolate temperature
    t_interp = t1 + t * (t2 - t1)
    
    return t_interp
end

function update_fixed_hsource_on_edge!(ps::KB.AbstractPhysicalSpace2D, ctr::KB.AM, 
    pos::Matrix{Float64}, area::Vector{Float64}, temperature::Vector{Float64}, hsource_1bc::Matrix{Float64}, ib::IBM2D)
    flags = ib.flags
    edges, nodes, pd_ids, xbis, nbis, xips, ip_nids, ip_bids = 
        ib.f2pd.edges, ib.f2pd.nodes, ib.f2pd.idpd, ib.f2pd.xbi, ib.f2pd.nbi, ib.f2pd.xip, ib.f2pd.idn, ib.f2pd.idb

    tem_ref = ib.rv.temperature
    l_ref = ib.rv.length
    k = ib.rv.k
    
    for i in axes(hsource_1bc, 2)
        hsource_1bc[i] = 0.0
    end

    for iter in axes(edges, 1)
        edge = edges[iter]
        x, y = xips[iter]
        node = nodes[iter]
        pd_id = pd_ids[iter]
        pf = [1, x, y, x * y]
        nx, ny = -1.0 .* xbis[iter]
        nids = ip_nids[iter]
        bids = ip_bids[iter]
        t1 = 0.5 * (temperature[pd_id[1][1]] + temperature[pd_id[1][2]]) / tem_ref 
        t2 = 0.5 * (temperature[pd_id[2][1]] + temperature[pd_id[2][2]]) / tem_ref

        hsource_1bc[node] += (100.0 * norm(edge[1] .- edge[2]) / area[node]) * (k * tem_ref / l_ref)
    end
end