using KitBase, OffsetArrays

function build_mesh(x0, x1, nx, y0, y1, ny, X, Y)
    mesh_x = zeros(size(X, 1)-1, size(Y, 1)-1)
    mesh_y = zeros(size(X, 1)-1, size(Y, 1)-1)
    dx = zeros(size(X, 1)-1, size(Y, 1)-1)
    dy = zeros(size(X, 1)-1, size(Y, 1)-1)
    vertices = zeros(size(X, 1)-1, size(Y, 1)-1, 4, 2)
    areas = zeros(size(X, 1)-1, size(Y, 1)-1, 4)
    n = Array{Vector{Float64}, 3}(undef, size(X, 1)-1, size(Y, 1)-1, 4)

    for i in 1:size(X, 1)-1, j in 1:size(Y, 1)-1
        mesh_x[i, j] = (X[i] + X[i+1]) / 2
        mesh_y[i, j] = (Y[j] + Y[j+1]) / 2
        dx[i, j] = X[i+1] - X[i]
        dy[i, j] = Y[j+1] - Y[j]
        n[i, j, 1] = [0.0, -1.0]
        n[i, j, 2] = [1.0, 0.0]
        n[i, j, 3] = [0.0, 1.0]
        n[i, j, 4] = [-1.0, 0.0]
        areas[i, j, 1] = dx[i, j]
        areas[i, j, 2] = dy[i, j]
        areas[i, j, 3] = dx[i, j]
        areas[i, j, 4] = dy[i, j]
        vertices[i, j, 1, 1] = X[i]
        vertices[i, j, 1, 2] = Y[j]
        vertices[i, j, 2, 1] = X[i+1]
        vertices[i, j, 2, 2] = Y[j]
        vertices[i, j, 3, 1] = X[i+1]
        vertices[i, j, 3, 2] = Y[j+1]
        vertices[i, j, 4, 1] = X[i]
        vertices[i, j, 4, 2] = Y[j+1]
    end

    ps = PSpace2D(x0, x1, nx, y0, y1, ny, 
                OffsetArray(mesh_x, 0:size(mesh_x, 1)-1, 0:size(mesh_x, 2)-1), 
                OffsetArray(mesh_y, 0:size(mesh_y, 1)-1, 0:size(mesh_y, 2)-1), 
                OffsetArray(dx, 0:size(dx, 1)-1, 0:size(dx, 2)-1), 
                OffsetArray(dy, 0:size(dy, 1)-1, 0:size(dy, 2)-1), 
                OffsetArray(vertices, 0:size(vertices, 1)-1, 0:size(vertices, 2)-1, 1:size(vertices, 3), 1:size(vertices, 4)), 
                OffsetArray(areas, 0:size(areas, 1)-1, 0:size(areas, 2)-1, 1:size(areas, 3)), 
                OffsetArray(n, 0:size(n, 1)-1, 0:size(n, 2)-1, 1:size(n, 3)))

    return ps
end