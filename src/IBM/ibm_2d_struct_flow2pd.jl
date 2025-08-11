mutable struct FluidToPD{T1, T2, T3, T4, T5}
    edges::T1   # Coordinates of the two endpoints of boundary segments
    nodes::T2   # Indices of PD points corresponding to the boundary
    idpd::T3    # Indices of PD points corresponding to segment endpoints (for solid solver)
    xbi::T4     # Coordinates of boundary intersection (bi) points
    nbi::T4     # Normal vectors at bi points
    xip::T4     # Mirror points of PD points about the boundary
    idn::T5     # Fluid grid indices corresponding to ip points (for fluid solver)
    idb::T5     # Solid grid indices corresponding to ip points (for fluid solver)

    function FluidToPD(ps::KB.AbstractPhysicalSpace2D, pos::Matrix{Float64}, nodes::Vector{Int}, boundarys::Vector{Vector{Vector{Float64}}}, bm2pd::Dict{Vector{Float64}, Vector{Int}}, flags::KB.AM)
        idpds = [Vector{Vector{Int}}() for iter = 1:size(boundarys, 1)]
        xbis = [Vector{Float64}() for iter = 1:size(boundarys, 1)]
        nbis = [Vector{Float64}() for iter = 1:size(boundarys, 1)]
        xips = [Vector{Float64}() for iter = 1:size(boundarys, 1)]
        ip_nids = [Vector{CartesianIndex}() for iter = 1:size(boundarys, 1)]
        ip_bids = [Vector{CartesianIndex}() for iter = 1:size(boundarys, 1)]
        
        for iter in axes(boundarys, 1)
            inside_point = pos[1:2, nodes[iter]]
            idpds[iter] = [bm2pd[boundarys[iter][1]], bm2pd[boundarys[iter][2]]]
            xbis[iter] = foot(inside_point, boundarys[iter][1], boundarys[iter][2])
            nbis[iter] = (xbis[iter] .- inside_point) / norm(xbis[iter] .- inside_point)
            xips[iter] = xbis[iter] .+ norm(xbis[iter] .- inside_point) * nbis[iter]
        end

        ip_nids, ip_bids = ip_connectivity(ps, xips, boundarys, flags)

        new{
            typeof(boundarys), 
            typeof(nodes),  
            typeof(idpds), 
            typeof(xbis),
            typeof(ip_bids),
        }(
            boundarys, # idg
            nodes,
            idpds,
            xbis, # xb
            nbis, # nb
            xips, # xi
            ip_nids, # idin
            ip_bids, # idib
        )

    end
end

function ip_connectivity(ps::KB.AbstractPhysicalSpace2D, xips::Vector{Vector{Float64}}, xbms::Vector{Vector{Vector{Float64}}}, flags::KB.AM)
    ip_nids = [Vector{CartesianIndex}() for iter = 1:size(xips, 1)]
    ip_bids = [Vector{CartesianIndex}() for iter = 1:size(xips, 1)]

    for iter in axes(xips, 1)
        xip = xips[iter]
        x, y = xips[iter]

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
                if are_points_on_the_same_side(xip, [x1, y1], xbms[iter])
                    push!(nids, CartesianIndex(idx + i, idy + j))
                else
                    push!(bids, CartesianIndex(idx + i, idy + j))
                end
            end
        end

        ip_nids[iter] = nids
        ip_bids[iter] = bids
    end

    return ip_nids, ip_bids
end