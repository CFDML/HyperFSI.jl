mutable struct IBM2D
    flags::KB.AM
    f2pd::FluidToPD
    pd2f::PDToFluid
    rv::ReferenceVariables

    function IBM2D(ps::KB.AbstractPhysicalSpace2D, ini_pos::Matrix{Float64}, ini_bc_edge::Dict{Int64, Vector{Vector{Vector{Float64}}}}, rv::ReferenceVariables)
        boundarys = Vector{Vector{Vector{Float64}}}()
        nodes = Vector{Int}()
        bm2pd = Dict{Vector{Float64}, Vector{Int}}()

        trans_para = rv.geo
        pos = compute_similarity_transform(ini_pos, trans_para)
        bc_edge = compute_similarity_transform(ini_bc_edge, trans_para)

        for (key, lines) in bc_edge
            for line in lines
                xbm1 = line[1][1:2]
                xbm2 = line[2][1:2]

                if haskey(bm2pd, xbm1)
                    push!(bm2pd[xbm1], key)
                else
                    bm2pd[xbm1] = []
                    push!(bm2pd[xbm1], key)
                end

                if haskey(bm2pd, xbm2)
                    push!(bm2pd[xbm2], key)
                else
                    bm2pd[xbm2] = []
                    push!(bm2pd[xbm2], key)
                end

                push!(boundarys, [xbm1, xbm2])
                push!(nodes, key)
            end
        end

        flags = build_flags(ps, boundarys)

        pd2flow = PDToFluid(ps, boundarys, bm2pd, flags)
        flow2pd = FluidToPD(ps, pos, nodes, boundarys, bm2pd, flags)

        new(flags, flow2pd, pd2flow, rv)
    end
end