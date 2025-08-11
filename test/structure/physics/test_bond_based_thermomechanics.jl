@testitem "BBTMMaterial_2D" begin
    using  HyperFSI: BBTMMaterial
    standard_thermomech_params = Dict{Symbol, Float64}(
        :k => 1.0,
        :cv => 1.0,
        :rho => 1.0,
        :δ => 2.0,
        :E => 2.0e11,
        :ε_c => 1e-2,
        :aph => 1.0e-6,
        :cv => 1.0,
        :rft => 273.0,
        :h => 25.0,
        :hσ => 5.73e-8,
        :hϵ => 0.9,
        :tem∞ => 0.0,
        :thick => 0.5
    )
    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5
    body = Peridynamics.Body(BBTMMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height) 

    @test body.mat isa BBTMMaterial
    @test body.point_params[1].rho == standard_thermomech_params[:rho]
    @test body.point_params[1].kc == standard_thermomech_params[:k]
    @test body.point_params[1].cv == standard_thermomech_params[:cv]
    @test body.point_params[1].δ == standard_thermomech_params[:δ]
    @test body.point_params[1].E == standard_thermomech_params[:E]
    @test body.point_params[1].εc == standard_thermomech_params[:ε_c]
    @test body.point_params[1].aph == standard_thermomech_params[:aph]
    @test body.point_params[1].rft == standard_thermomech_params[:rft]
    @test body.point_params[1].h == standard_thermomech_params[:h]
    @test body.point_params[1].hσ == standard_thermomech_params[:hσ]
    @test body.point_params[1].hϵ == standard_thermomech_params[:hϵ]
    @test body.point_params[1].tem∞ == standard_thermomech_params[:tem∞]
    @test body.point_params[1].nu == 1/3
    @test body.point_params[1].kp == 6*standard_thermomech_params[:k] / (π * height * standard_thermomech_params[:δ]^3)
    @test body.point_params[1].bc == 9*standard_thermomech_params[:E] / (π * height * standard_thermomech_params[:δ]^3)
end

@testitem "BBTMaterial_3D" begin
    using  HyperFSI: BBTMMaterial
    standard_thermomech_params = Dict{Symbol, Float64}(
        :k => 1.0,
        :cv => 1.0,
        :rho => 1.0,
        :δ => 2.0,
        :E => 2.0e11,
        :ε_c => 1e-2,
        :aph => 1.0e-6,
        :cv => 1.0,
        :rft => 273.0,
        :h => 25.0,
        :hσ => 5.73e-8,
        :hϵ => 0.9,
        :tem∞ => 0.0,
    )
    pos  = [ -0.25   0.25  -0.25   0.25  -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25   0.25  -0.25  -0.25   0.25  0.25
             -0.25  -0.25  -0.25  -0.25   0.25   0.25   0.25  0.25]
    vol = [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125]

    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0) 

    @test body.mat isa BBTMaterial
    @test body.point_params[1].rho == standard_thermomech_params[:rho]
    @test body.point_params[1].kc == standard_thermomech_params[:k]
    @test body.point_params[1].cv == standard_thermomech_params[:cv]
    @test body.point_params[1].δ == standard_thermomech_params[:δ]
    @test body.point_params[1].E == standard_thermomech_params[:E]
    @test body.point_params[1].εc == standard_thermomech_params[:ε_c]
    @test body.point_params[1].aph == standard_thermomech_params[:aph]
    @test body.point_params[1].rft == standard_thermomech_params[:rft]
    @test body.point_params[1].h == standard_thermomech_params[:h]
    @test body.point_params[1].hσ == standard_thermomech_params[:hσ]
    @test body.point_params[1].hϵ == standard_thermomech_params[:hϵ]
    @test body.point_params[1].tem∞ == standard_thermomech_params[:tem∞]
    @test body.point_params[1].nu == 1/4
    @test body.point_params[1].kp == 6*standard_thermomech_params[:k] / (π * standard_thermomech_params[:δ]^4)
    @test body.point_params[1].bc == 12*standard_thermomech_params[:E] / (π * standard_thermomech_params[:δ]^4)
end
