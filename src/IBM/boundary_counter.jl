
mutable struct BoundaryCounter
    counts::Dict{Tuple{Vararg{Int}}, Int}                    # Key without direction → count
    records::Dict{Tuple{Vararg{Int}}, Vector{Tuple{Vector{Int},Int}}}  # Key without direction → [(face_origional, pd_id)]
end

function BoundaryCounter()
    BoundaryCounter(Dict{Tuple{Vararg{Int}}, Int}(),
                    Dict{Tuple{Vararg{Int}}, Vector{Tuple{Vector{Int},Int}}}())
end

function add_element!(bc::BoundaryCounter, element::Vector{Int}, celltype::String ,pd_id::Int)
    for face in extract_faces(element, celltype)  # Assuming 2D triangles; modify as needed
        key = sort(face) |> Tuple
        bc.counts[key] = get(bc.counts, key, 0) + 1
        recs = get!(bc.records, key, Vector{Tuple{Vector{Int},Int}}())
        push!(recs, (face, pd_id))  # Keep original direction
    end
end

function remove_element!(bc::BoundaryCounter, element::Vector{Int}, celltype::String, pd_id::Int)
    for face in extract_faces(element, celltype)
        key = sort(face) |> Tuple
        if !haskey(bc.counts, key)
            error("remove_element!: face $face not found")
        end

        bc.counts[key] -= 1
        recs = bc.records[key]

        for i in eachindex(recs)
            if recs[i][1] == face && recs[i][2] == pd_id
                recs[i] = recs[end]  
                pop!(recs)
                break
            end
        end

        if bc.counts[key] == 0
            delete!(bc.counts, key)
            delete!(bc.records, key)
        end
    end
end

function get_boundary(bc::BoundaryCounter)
    result = Vector{Tuple{Vector{Int},Int}}() 
    for (key,cnt) in bc.counts
        if cnt == 1
            append!(result, bc.records[key])
        end
    end
    return result
end

function extract_faces(element::Vector{Int}, celltype::String)
    n = length(element)
    celltype_upper = uppercase(celltype)

    if celltype_upper == "CPS3" && n == 3
        return [[element[1],element[2]],
                [element[2],element[3]],
                [element[3],element[1]]]

    elseif celltype_upper == "CPS4" && n == 4
        return [[element[1],element[2]],
                [element[2],element[3]],
                [element[3],element[4]],
                [element[4],element[1]]]

    elseif celltype_upper == "C3D4" && n == 4
        return [[element[1],element[2],element[3]],
                [element[1],element[4],element[2]],
                [element[2],element[4],element[3]],
                [element[3],element[4],element[1]]]

    elseif celltype_upper == "C3D8" && n == 8
        return [[element[1],element[2],element[3],element[4]],
                [element[5],element[6],element[7],element[8]],
                [element[1],element[5],element[6],element[2]],
                [element[2],element[6],element[7],element[3]],
                [element[3],element[7],element[8],element[4]],
                [element[4],element[8],element[5],element[1]]]

    # 3D（Prism）
    elseif (celltype_upper == "C3D6" || celltype_upper == "PRISM") && n == 6
        return [
            [element[1], element[2], element[3]],  # Bottom
            [element[4], element[5], element[6]],  # Top
            [element[1], element[2], element[5], element[4]],  # Side 1
            [element[2], element[3], element[6], element[5]],  # Side 2
            [element[3], element[1], element[4], element[6]]   # Side 3
        ]

    else
        error("Unsupported cell type or element length mismatch: $celltype with $n nodes")
    end
end
