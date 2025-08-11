function build_flags(ps::KB.AbstractPhysicalSpace2D, boundarys::Vector{Vector{Vector{Float64}}})
    flags = ones(Int, axes(ps.x))
    for i in axes(flags, 1), j in axes(flags, 2)
        point = [ps.x[i,j], ps.y[i,j]]
        if is_point_in_polygon(point, boundarys)
            flags[i, j] = 0
        end
    end
    flags[0, :] .= -1       # left
    flags[ps.nx+1, :] .= -1 # right
    flags[:, 0] .= -1       # down
    flags[:, ps.ny+1] .= -1 # up

    KB.ghost_flag!(ps, flags)
    
    return flags
end

function compute_similarity_transform!(data::Dict{Int, Vector{Vector{Vector{Float64}}}}, trans_para::Vector{Float64})
    length(trans_para) == 4 || error("变换参数需要 [s, θ, tx, ty]")
    s, θ, tx, ty = trans_para
    
    R = [cos(θ) -sin(θ) 0;
         sin(θ)  cos(θ) 0;
         0       0      1]
    
    for (_, edges) in data
        for edge in edges
            for point in edge
                v = [point[1], point[2], point[3]]
                v_transformed = s * R * v + [tx, ty, 0.0]
                point[1], point[2], point[3] = v_transformed[1], v_transformed[2], v_transformed[3]
            end
        end
    end
    return data
end

function compute_similarity_transform(data::Dict{Int, Vector{Vector{Vector{Float64}}}}, trans_para::Vector{Float64})
    new_data = deepcopy(data)
    compute_similarity_transform!(new_data, trans_para)
    return new_data
end

function compute_similarity_transform(points::Matrix{Float64}, trans_para::Vector{Float64})
    length(trans_para) == 4 || error("变换参数需要 [s, θ, tx, ty]")
    s, θ, tx, ty = trans_para

    R = [cos(θ) -sin(θ) 0;
         sin(θ)  cos(θ) 0;
         0       0      1]
    
    sR = s * R
    transformed = sR * points .+ [tx, ty, 0]

    return transformed
end


function is_point_in_polygon(point::Vector{Float64}, polygon::Vector{Vector{Vector{Float64}}})
    x, y = point

    n = length(polygon)
    inside = false

    for i in 1:n
        x1, y1 = polygon[i][1]
        x2, y2 = polygon[i][2]
        
        if is_point_on_line_segment(point, polygon[i][1], polygon[i][2])
            return true
        end

        if ((y1 > y) != (y2 > y)) && (x < (x2 - x1) * (y - y1) / (y2 - y1) + x1)
            inside = !inside
        end
    end
    
    return inside
end

function is_point_on_line_segment(p::Vector{Float64}, p1::Vector{Float64}, p2::Vector{Float64})::Bool
    x, y = p
    x1, y1 = p1
    x2, y2 = p2
    
    if p1 == p2
        return p == p1
    end


    cross_product = (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)
    
    if cross_product != 0
        return false
    end

    if min(x1, x2) <= x <= max(x1, x2) && min(y1, y2) <= y <= max(y1, y2)
        return true
    else
        return false
    end
end

function are_points_on_the_same_side(A::Vector{Float64}, B::Vector{Float64}, l::Vector{Vector{Float64}})::Bool
    P = l[1]
    Q = l[2]

    PQ = Q - P
    PA = A - P
    PB = B - P

    cross_p = PQ[1] * PA[2] - PQ[2] * PA[1]
    cross_q = PQ[1] * PB[2] - PQ[2] * PB[1]

    if cross_p * cross_q > 0
        return true
    else
        return false
    end
    return true
end

function foot(p::Vector{Float64}, p1::Vector{Float64}, p2::Vector{Float64})

    A = p2 .- p1
    B = p .- p1

    t = dot(A, B) / dot(A, A)

    foot = p1 .+ t .* A
    
    return foot
end

struct LineSegment
    p1::Vector{Float64}
    p2::Vector{Float64}
end

function orientation(p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64})
    val = (p2[2] - p1[2]) * (p3[1] - p2[1]) - (p2[1] - p1[1]) * (p3[2] - p2[2])
    if val == 0
        return 0  # collinear
    elseif val > 0
        return 1  # clockwise
    else
        return 2  # counterclockwise
    end
end

function on_segment(p::Vector{Float64}, q::Vector{Float64}, r::Vector{Float64})
    if q[1] <= max(p[1], r[1]) && q[1] >= min(p[1], r[1]) && q[2] <= max(p[2], r[2]) && q[2] >= min(p[2], r[2])
        return true
    end
    return false
end

function do_intersect(s1::LineSegment, s2::LineSegment)::Bool
    p1, q1 = s1.p1, s1.p2
    p2, q2 = s2.p1, s2.p2

    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    if o1 != o2 && o3 != o4
        return true
    end

    if o1 == 0 && on_segment(p1, p2, q1)
        return true
    end

    if o2 == 0 && on_segment(p1, q2, q1)
        return true
    end

    if o3 == 0 && on_segment(p2, p1, q2)
        return true
    end

    if o4 == 0 && on_segment(p2, q1, q2)
        return true
    end

    return false
end

function do_intersect(poly_vertices::Vector{Vector{Vector{Float64}}}, inside_point::Vector{Float64}, outside_point::Vector{Float64})
    
    s2 = LineSegment(inside_point, outside_point)

    for i in 1:size(poly_vertices, 1)
        s1 = LineSegment(poly_vertices[i][1], poly_vertices[i][2])
        if do_intersect(s1, s2)
            return i
        end
    end

    error("$(inside_point) 和 $(outside_point) 在多边形同一侧")
end