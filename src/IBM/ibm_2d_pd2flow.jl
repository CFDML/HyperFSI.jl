function update_ghost!(ctr::KB.AM, ps::KitBase.AbstractPhysicalSpace2D, gas::Gas, ib::IBM2D, temperature::Vector{Float64})
    ghost_ids, xbms, pd_ids, xbis, nbis, xips, ip_nids, ip_bids, ghost_ctr = 
        ib.pd2f.idgc, ib.pd2f.xbm, ib.pd2f.idpd, ib.pd2f.xbi, ib.pd2f.nbi, ib.pd2f.xip, ib.pd2f.idn, ib.pd2f.idb, ib.pd2f.ctr

        flags = ib.flags
    tem_ref = ib.rv.temperature

    @threads for iter in axes(ghost_ids, 1)
        for (direc, nid) in ip_nids[iter]
            bid = ip_bids[iter][direc]
            xbi = xbis[iter][direc]
            xf, yf = xips[iter][direc]                
            pos = [1, xf, yf, xf * yf]
            xbm = xbms[iter][direc]
            pd_id = pd_ids[iter][direc]
            t1 = 0.5 * (temperature[pd_id[1][1]] + temperature[pd_id[1][2]]) / tem_ref
            t2 = 0.5 * (temperature[pd_id[2][1]] + temperature[pd_id[2][2]]) / tem_ref
            if size(nid, 1) > 1 && size(nid, 1) + size(bid, 1) == 4
                # U
                w1 = [ctr[idx].prim[2] for idx in nid]
                w2 = zeros(size(bid, 1))
                w = [w1; w2]

                C = bilinear_coeffs(ps, xbm, nid, bid, w)
                U1 = C' * pos

                # V
                w1 = [ctr[idx].prim[3] for idx in nid]
                w2 = zeros(size(bid, 1))
                w = [w1; w2]

                C = bilinear_coeffs(ps, xbm, nid, bid, w)
                V1 = C' * pos

                # T
                w1 = [1 / ctr[idx].prim[4] for idx in nid]
                w2 = ones(size(bid, 1))
                for i in axes(bid, 1)
                    xb, yb = foot([ps.x[bid[i]], ps.y[bid[i]]], xbm[1], xbm[2])
                    w2[i] = temperature_interpolation([xb, yb], xbm[1], t1, xbm[2], t2)
                end
                w = [w1; w2]

                C = bilinear_coeffs(ps, xbm, nid, bid, w)
                T1 = C' * pos

                # p
                w1 = [0.5 * ctr[idx].prim[1] / ctr[idx].prim[4] for idx in nid]
                w = w1

                C = bilinear_coeffs(ps, xbm, nid, bid, w)
                P1 = C' * pos

                ρ1 = 2 * P1 / T1
            elseif size(nid, 1) >= 1
                idx1 = nid[1]
                ρ1, U1, V1, λ1 = ctr[idx1].prim
                T1 = 1.0 / λ1
                P1 = 0.5 * ρ1 / λ1
            else
                error("No fluid nodes found around the IP node!")
            end

            tp = temperature_interpolation(xbi, xbm[1], t1, xbm[2], t2)
            T0 = tp * 2 - T1
            ρ0 = ρ1 * T1 / T0
            U0 = -U1
            V0 = -V1

            ghost_ctr[iter][direc].prim .= [ρ0, U0, V0, 1 / T0]
            ghost_ctr[iter][direc].w .= prim_conserve(ghost_ctr[iter][direc].prim, gas.γ)
        end
    end
    return nothing
end