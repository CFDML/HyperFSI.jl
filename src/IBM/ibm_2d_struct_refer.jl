struct ReferenceVariables
    length::Float64
    velocity::Float64
    time::Float64
    density::Float64
    temperature::Float64
    pressure::Float64
    R::Float64
    μ::Float64
    k::Float64
    geo::Vector{Float64}

    function ReferenceVariables(len, vel, time, density, temp, press, R, μ, k, geo)
        len > 0 || error("len must be positive")
        density > 0 || error("density must be positive")
        temp > 0 || error("temp must be positive")
        press > 0 || error("press must be positive")
        new(len, vel, time, density, temp, press, R, μ, k, geo)
    end
end

function ReferenceVariables(;
    length = -1.0,
    velocity = -1.0,
    time = -1.0,
    density = -1.0,
    temperature = -1.0,
    pressure = -1.0,
    R = -1.0,
    m = -1.0,
    μ = -1.0,
    k = -1.0,
    geo = nothing,
)
    if velocity < 0.0
        if R > 0
            velocity = (2 * R * temperature) ^ 0.5
        elseif m > 0
            velocity = (2 * 8.314 * temperature / m) ^ 0.5
        else
            error("Gas constant or molecular mass must be specified to compute reference velocity")
        end
    end

    if R < 0.0 && m > 0.0
        R = 8.31 / m
    elseif R < 0.0 && m < 0.0
        m = 29 * 0.001
        R = 8.31 / m
    end

    if μ < 0.0
        μ = density * velocity * length
    end

    if k < 0.0
        k = μ * 3.5 * R
    end

    if time < 0
        time = length / velocity
    end

    if isnothing(geo)
        geo = [1.0, 0.0, 0.0, 0.0]
    end

    ReferenceVariables(length, velocity, time, density, temperature, pressure, R, μ, k, geo)
end