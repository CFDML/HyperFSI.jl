@testitem "BBTMaterial_2D" begin
    using  HyperFSI: BBTMaterial
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

    @test body.mat isa BBTMaterial
    @test body.point_params[1].rho == standard_thermal_params[:rho]
    @test body.point_params[1].kc == standard_thermal_params[:k]
    @test body.point_params[1].cv == standard_thermal_params[:cv]
    @test body.point_params[1].δ == standard_thermal_params[:δ]
    @test body.point_params[1].kp == 6*standard_thermal_params[:k] / (π * height * standard_thermal_params[:δ]^3)
end

@testitem "BBTMaterial_3D" begin
    using  HyperFSI: BBTMaterial
    standard_thermal_params = Dict{Symbol, Float64}(
        :k => 1.0,
        :cv => 1.0,
        :rho => 1.0,
        :δ => 2.0
    )
    pos  = [ -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
             -0.25  -0.25  -0.25  -0.25   0.25   0.25   0.25  0.25]
    vol = [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125]

    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0) 

    @test body.mat isa BBTMaterial
    @test body.point_params[1].rho == standard_thermal_params[:rho]
    @test body.point_params[1].kc == standard_thermal_params[:k]
    @test body.point_params[1].cv == standard_thermal_params[:cv]
    @test body.point_params[1].δ == standard_thermal_params[:δ]
    @test body.point_params[1].kp == 6*standard_thermal_params[:k] / (π * standard_thermal_params[:δ]^4)
end
