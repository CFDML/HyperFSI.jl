function export_fluid(ps::KitBase.AbstractPhysicalSpace2D, ctr::KB.AM, options::Peridynamics.AbstractJobOptions, tt::Int, Δt::Float64)
    file_path = options.flow
    if !isdir(file_path)
        mkdir(file_path)
    end
    if mod(tt, options.freq) == 0
        write_fluid_vtk(ps, ctr, file_path, tt, Δt)
    end
    return nothing
end

function write_fluid_vtk(ps, ctr, file_path, tt, dtf)
    # x = ps.x[:, 1]
    # y = ps.y[1, :]
    x = [ps.x[i, 1] - ps.dx[i, 1]/2 for i in axes(ps.x, 1)]
    push!(x, ps.x[end, 1]+ps.dx[end, 1]/2)
    y = [ps.y[1, i] - ps.dy[1, i]/2 for i in axes(ps.x, 2)]
    push!(y, ps.y[1, end]+ps.dy[1, end]/2)

    velocity = Array{Float64}(undef, 2, size(ctr, 1), size(ctr, 2))
    for i in axes(ctr, 1), j in axes(ctr, 2)
        velocity[1, i+1, j+1] = ctr[i,j].prim[2]
        velocity[2, i+1, j+1] = ctr[i,j].prim[3]
    end

    saved_files = paraview_collection(joinpath(file_path, "full_simulation"); append = true) do pvd
        vtk_grid(joinpath(file_path, "fluid_output_$(tt)"), x, y) do vtk
            vtk["Density"] = [ctr[i,j].prim[1] for i in axes(ps.x,1), j in axes(ps.x, 2)]
            vtk["Velocity"] = velocity
            vtk["Temperature"] = [1.0 / ctr[i,j].prim[end] for i in axes(ps.x,1), j in axes(ps.x, 2)]
            vtk["Pressure"] = [0.5 * ctr[i,j].prim[1] / ctr[i,j].prim[end] for i in axes(ps.x,1), j in axes(ps.x, 2)]
            pvd[tt*dtf] = vtk
        end  
    end
end

function output_to_main!(local_vars::Dict{Symbol, Any}, options::Peridynamics.AbstractJobOptions, tt::Int, output, output_vars)
    if !isnothing(output_vars) && mod(tt, options.freq) == 0
        A = []
        for a in output_vars
            if haskey(local_vars, Symbol(a))
                push!(A, deepcopy(local_vars[Symbol(a)]))
            else
                error("变量 '$b' 不存在")
            end
        end
        push!(output, A)
    end
end