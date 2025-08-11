mutable struct PDToFluid{TC,TF,T1,T2,T3,T4}
    idgc::TC   # Grid index of GC (fluid solver)  
    xbm::TF    # Coordinates of BM points at segment ends  
    idpd::T1   # PD point indices for segments (solid solver)  
    xbi::T2    # BI point coordinates  
    nbi::T2    # Normal vectors at BI points  
    xip::T2    # IP point coordinates  
    idn::T3    # Fluid grid indices for IP points (fluid solver)  
    idb::T3    # Solid grid indices for IP points (fluid solver)  
    ctr::T4    # Pseudo-flow variables stored in GC  

    function PDToFluid(ps::KB.AbstractPhysicalSpace2D, boundarys::Vector{Vector{Vector{Float64}}}, bm2pd::Dict{Vector{Float64}, Vector{Int}}, flags::KB.AM)
        direcs = [CartesianIndex(-1, 0), CartesianIndex(1, 0), CartesianIndex(0, -1), CartesianIndex(0, 1)]
    
        ghost_ids = findall(flags .== -2)
        xbms = [Dict{CartesianIndex, Vector{Vector{Float64}}}() for iter = 1:size(ghost_ids, 1)]
        idpds = [Dict{CartesianIndex, Vector{Vector{Int}}}() for iter = 1:size(ghost_ids, 1)]
        xbis = [Dict{CartesianIndex, Vector{Float64}}() for iter = 1:size(ghost_ids, 1)]
        nbis = [Dict{CartesianIndex, Vector{Float64}}() for iter = 1:size(ghost_ids, 1)]
        xips = [Dict{CartesianIndex, Vector{Float64}}() for iter = 1:size(ghost_ids, 1)]
        ip_nids = [Dict{CartesianIndex, Vector{CartesianIndex}}() for iter = 1:size(ghost_ids, 1)]
        ip_bids = [Dict{CartesianIndex, Vector{CartesianIndex}}() for iter = 1:size(ghost_ids, 1)]
    
        for iter in axes(xbis, 1)
            idx = ghost_ids[iter]
    
            inside_point = [ps.x[idx], ps.y[idx]]
            for direc in direcs
                target_cell = idx + direc
                if flags[target_cell] == 1
                    outside_point = [ps.x[target_cell], ps.y[target_cell]]
                    idx1 = do_intersect(boundarys, inside_point, outside_point)
                    xbms[iter][direc] = boundarys[idx1]
                    idpds[iter][direc] = [bm2pd[xbms[iter][direc][1]],bm2pd[xbms[iter][direc][2]]]
                    xbis[iter][direc] = foot(inside_point, boundarys[idx1][1], boundarys[idx1][2])
                    nbis[iter][direc] = (xbis[iter][direc] .- inside_point) / norm(xbis[iter][direc] .- inside_point)
                    xips[iter][direc] = xbis[iter][direc] .+ norm(xbis[iter][direc] .- inside_point) * nbis[iter][direc]
                end
            end        
        end
    
        ip_nids, ip_bids = ip_connectivity(ps, xips, xbms, flags)
    
        ghost_ctr = [Dict{CartesianIndex, ControlVolume}() for iter = 1:size(ghost_ids, 1)]
    
        for iter in axes(ghost_ctr, 1)
            for (direc, xbi) in xbis[iter]
                ghost_ctr[iter][direc] = ControlVolume(zeros(Float64, 4), zeros(Float64, 4), 2)
            end
        end

        new{
            typeof(ghost_ids), 
            typeof(xbms),  
            typeof(idpds), 
            typeof(xbis),
            typeof(ip_bids),
            typeof(ghost_ctr),
        }(
            ghost_ids, # idg
            xbms,
            idpds,
            xbis, # xb
            nbis, # nb
            xips, # xi
            ip_nids, # idin
            ip_bids, # idib
            ghost_ctr,
        )
    end
end

function ip_connectivity(ps::KB.AbstractPhysicalSpace2D, xips::Vector{Dict{CartesianIndex, Vector{Float64}}}, xbms::Vector{Dict{CartesianIndex, Vector{Vector{Float64}}}}, flags::KB.AM)
    ip_nids = [Dict{CartesianIndex, Vector{CartesianIndex}}() for iter = 1:size(xips, 1)]
    ip_bids = [Dict{CartesianIndex, Vector{CartesianIndex}}() for iter = 1:size(xips, 1)]

    for iter in eachindex(xips)
        for (direc, xip) in xips[iter]
            x, y = xip

            # id of the cell where 
            dxs = abs.(x .- ps.x[:, 1])
            dys = abs.(y .- ps.y[1, :])
            cidx = argmin(dxs)
            cidy = argmin(dys)

            # id of the center face intersection of the interpolation stencil
            idx = begin
                if x > ps.x[cidx, 1]
                    cidx + 1
                else
                    cidx
                end
            end
            idy = begin
                if y > ps.y[1, cidy]
                    cidy + 1
                else
                    cidy
                end
            end

            nids = CartesianIndex[]
            bids = CartesianIndex[]

            for i in -1:0, j in -1:0
                x1 = ps.x[idx + i, idy + j]
                y1 = ps.y[idx + i, idy + j]
                if flags[idx + i, idy + j] != -1
                    if are_points_on_the_same_side(xip, [x1, y1], xbms[iter][direc])
                        push!(nids, CartesianIndex(idx + i, idy + j))
                    else
                        push!(bids, CartesianIndex(idx + i, idy + j))
                    end
                end
            end

            ip_nids[iter][direc] = nids
            ip_bids[iter][direc] = bids
        end
    end

    return ip_nids, ip_bids
end


