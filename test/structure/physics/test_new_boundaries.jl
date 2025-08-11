@testitem "temperature_ic!_single_dim_ics" begin
    using  HyperFSI: temperature_ic!
    standard_thermal_params = Dict{Symbol, Float64}(
        :k => 1.0,
        :cv => 1.0,
        :rho => 1.0,
        :δ => 2.0
    )
    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5
    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height) 

    temperature_ic!(body, :all_points, 100.0)

    @test body.posdep_single_dim_ics == []

    @test body.single_dim_ics[1].dim == 1
    @test body.single_dim_ics[1].field == :temperature
    @test body.single_dim_ics[1].point_set == :all_points
    @test body.single_dim_ics[1].value == 100.0
end

@testitem "temperature_ic!_posdep_single_dim_ics" begin
    using  HyperFSI: temperature_ic!
    standard_thermal_params = Dict{Symbol, Float64}(
        :k => 1.0,
        :cv => 1.0,
        :rho => 1.0,
        :δ => 2.0
    )
    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5
    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height) 

    temperature_ic!(p -> p[1]*10+p[2]*10, body, :all_points)

    @test body.single_dim_ics == []
    @test body.posdep_single_dim_ics[1].dim == 1
    @test body.posdep_single_dim_ics[1].field == :temperature
    @test body.posdep_single_dim_ics[1].point_set == :all_points
    @test body.posdep_single_dim_ics[1].fun(pos[:,1]) == -5.0
end

@testitem "temperature_bc!_single_dim_bcs" begin
    using  HyperFSI: temperature_bc!
    standard_thermal_params = Dict{Symbol, Float64}(
        :k => 1.0,
        :cv => 1.0,
        :rho => 1.0,
        :δ => 2.0
    )
    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5
    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height) 

    temperature_bc!(t -> 99t, body, :all_points)

    @test body.posdep_single_dim_bcs == []
    @test body.single_dim_bcs[1].dim == 1
    @test body.single_dim_bcs[1].field == :temperature
    @test body.single_dim_bcs[1].point_set == :all_points
    @test body.single_dim_bcs[1].fun(2) == 198.0
end

@testitem "temperature_bc!_poedep_single_dim_bcs" begin
    using  HyperFSI: temperature_bc!
    standard_thermal_params = Dict{Symbol, Float64}(
        :k => 1.0,
        :cv => 1.0,
        :rho => 1.0,
        :δ => 2.0
    )
    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5
    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height) 

    temperature_bc!((p,t) -> p[1]*10*t+p[2]*1*t, body, :all_points)

    @test body.single_dim_bcs == []
    @test body.posdep_single_dim_bcs[1].dim == 1
    @test body.posdep_single_dim_bcs[1].field == :temperature
    @test body.posdep_single_dim_bcs[1].point_set == :all_points
    @test body.posdep_single_dim_bcs[1].fun(pos[:,2],2) == 4.5
end

@testitem "temperature_databc!" begin
    using  HyperFSI: temperature_databc!
    standard_thermal_params = Dict{Symbol, Float64}(
        :k => 1.0,
        :cv => 1.0,
        :rho => 1.0,
        :δ => 2.0
    )
    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5

    data = pos[1:1,:] .* 10 .+ 5 

    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height) 

    temperature_databc!(body, data, :all_points)

    @test body.single_dim_bcs == []
    @test body.posdep_single_dim_bcs == []

    @test body.data_bcs[1].dims == [1]
    @test body.data_bcs[1].field == :temperature
    @test body.data_bcs[1].point_set == :all_points
    @test body.data_bcs[1].data[1,3] == 2.5
end

@testitem "hsource_bc!_single_dim_bcs" begin
    using  HyperFSI: hsource_bc!
    standard_thermal_params = Dict{Symbol, Float64}(
        :k => 1.0,
        :cv => 1.0,
        :rho => 1.0,
        :δ => 2.0
    )
    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5
    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height) 

    hsource_bc!(t -> 1.1t, body, :all_points)

    @test body.posdep_single_dim_bcs == []
    @test body.single_dim_bcs[1].dim == 1
    @test body.single_dim_bcs[1].field == :hsource
    @test body.single_dim_bcs[1].point_set == :all_points
    @test body.single_dim_bcs[1].fun(2) == 2.2
end

@testitem "hsource_bc!_poedep_single_dim_bcs" begin
    using  HyperFSI: hsource_bc!
    standard_thermal_params = Dict{Symbol, Float64}(
        :k => 1.0,
        :cv => 1.0,
        :rho => 1.0,
        :δ => 2.0
    )
    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5
    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height) 

    hsource_bc!((p,t) -> p[1]*10*t+p[2]*10*t, body, :all_points)

    @test body.single_dim_bcs == []
    @test body.posdep_single_dim_bcs[1].dim == 1
    @test body.posdep_single_dim_bcs[1].field == :hsource
    @test body.posdep_single_dim_bcs[1].point_set == :all_points
    @test body.posdep_single_dim_bcs[1].fun(pos[:,2],2) == 0.0
end

@testitem "hsource_databc!" begin
    using  HyperFSI: temperature_databc!
    standard_thermal_params = Dict{Symbol, Float64}(
        :k => 1.0,
        :cv => 1.0,
        :rho => 1.0,
        :δ => 2.0
    )
    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5

    data = pos[1:1,:] .* 10 .+ 10 

    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height) 

    hsource_databc!(body, data, :all_points)

    @test body.single_dim_bcs == []
    @test body.posdep_single_dim_bcs == []

    @test body.data_bcs[1].dims == [1]
    @test body.data_bcs[1].field == :hsource
    @test body.data_bcs[1].point_set == :all_points
    @test body.data_bcs[1].data[1,3] == 7.5
end

@testitem "find_sec_bcs_points" begin
    using HyperFSI: find_sec_bcs_points
    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5
    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height) 

    Peridynamics.point_set!(p -> p[1]< 0.0, body, :radiation)
    Peridynamics.point_set!(p -> p[1]> 0.0, body, :convention)
    temperature_bc!(t -> 99t, body, :radiation)

    vv = Thermstep(steps=1) # Time solver setup for thermal diffusion simulation in structure
    job = Peridynamics.Job(body, vv)
    ret = FSI_submit(job, "T"; quiet=true)

    conv, radi = find_sec_bcs_points(ret)
   
    @test typeof(conv) == Vector{Vector{Int64}}
    @test typeof(radi) == Vector{Vector{Int64}}
    @test length(conv) == 1
    @test length(radi) == 1
    @test length(conv[1]) == 2
    @test length(radi[1]) == 2
    @test body.point_sets[:convention] == [2, 4]
    @test body.point_sets[:radiation] == [1, 3]
    @test conv[1] == [2, 4]
    @test radi[1] == [1, 3]
end

@testitem "second_bcs!" begin
    using HyperFSI: second_bcs!
    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5
    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=1e-8, hϵ=1.0, tem∞=0.0, thick=height) 

    Peridynamics.point_set!(p -> p[1]< 0.0, body, :radiation)
    Peridynamics.point_set!(p -> p[1]> 0.0, body, :convention)
    temperature_ic!(body, :all_points, 1000.0)

    vv = Thermstep(steps=1, stepsize = 1e-6) # Time solver setup for thermal diffusion simulation in structure
    job = Peridynamics.Job(body, vv)
    ret = FSI_submit(job, "T"; quiet=true)
   
    @test ret.chunks[1].storage.hsource[1, 1] == 1e-8*1*(0-1000^4)/0.5
    @test ret.chunks[1].storage.hsource[1, 2] == 25*(0-1000)/0.5
    @test ret.chunks[1].storage.hsource[1, 3] == 1e-8*1*(0-1000^4)/0.5
    @test ret.chunks[1].storage.hsource[1, 4] == 25*(0-1000)/0.5
    @test ret.chunks[1].storage.pflux[1, 1] == 0
    @test ret.chunks[1].storage.pflux[1, 2] == 0
    @test ret.chunks[1].storage.pflux[1, 3] == 0
    @test ret.chunks[1].storage.pflux[1, 4] == 0
    @test ret.chunks[1].storage.temperature[1, 1] == 999.98
    @test ret.chunks[1].storage.temperature[1, 2] == 999.95
    @test ret.chunks[1].storage.temperature[1, 3] == 999.98
    @test ret.chunks[1].storage.temperature[1, 4] == 999.95
end




