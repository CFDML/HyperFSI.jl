@testitem "BBTTMaterial_const_homogeneous" begin
    using  HyperFSI: BBTTMaterial

    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5
    body = Body(BBTTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height) 

    @test body.mat isa BBTTMaterial
    @test body.point_params[1].kp_T == HyperFSI.const_kp
    @test body.point_params[1].cv_T == HyperFSI.const_cv
    @test body.point_params[1].aph_T == HyperFSI.const_aph
    @test HyperFSI.const_kp(1.0, 2.0, 300.0) == 2.0
    @test HyperFSI.const_kp(1.0, 2.0, 1000.0) == HyperFSI.const_kp(1.0, 2.0, 300.0)
    @test HyperFSI.const_cv(1.0, 300.0) == 1.0
    @test HyperFSI.const_cv(1.0, 1000.0) == HyperFSI.const_cv(1.0, 300.0)
    @test HyperFSI.const_aph(1.0e-6, 300.0) == 1.0e-6
    @test HyperFSI.const_aph(1.0e-6, 1000.0) == HyperFSI.const_aph(1.0e-6, 300.0)
end

@testitem "BBTTMaterial_k(T)_homogeneous" begin
    using  HyperFSI: BBTTMaterial

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
    
    function k_linear(kc, kp, T) 
        kc_act = (1.0 + (T - 300)/300) # linear increase of k with T
        kp_act = kc_act * kp / kc
        return kp_act 
    end

    body = Body(BBTTMaterial(), pos, vol)
    material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, kp_T=k_linear)

    @test body.point_params[1].kp_T == k_linear
    @test body.point_params[1].cv_T == HyperFSI.const_cv
    @test body.point_params[1].aph_T == HyperFSI.const_aph
    @test body.point_params[1].kc == standard_thermal_params[:k]
    @test body.point_params[1].kp == 6*standard_thermal_params[:k] / (π * standard_thermal_params[:δ]^4)
    @test k_linear(1.0, body.point_params[1].kp, 300.0) == body.point_params[1].kp
    @test k_linear(1.0, body.point_params[1].kp, 600.0) == 2*body.point_params[1].kp
end

@testitem "BBTTMaterial_k(T)_cv(T)_homogeneous" begin
    using  HyperFSI: BBTTMaterial
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
    
    function k_linear(kc, kp, T) 
        kc_act = (1.0 + (T - 300)/300) # linear increase of k with T
        kp_act = kc_act * kp / kc
        return kp_act 
    end
    
    cv_linear(cv, T) = (1.0 +  (T - 300)/300) # linear increase of cv with T
    body = Body(BBTTMaterial(), pos, vol)
    material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, kp_T=k_linear, cv_T=cv_linear)

    @test body.point_params[1].kp_T == k_linear
    @test body.point_params[1].cv_T == cv_linear
    @test body.point_params[1].aph_T == HyperFSI.const_aph
    @test body.point_params[1].kc == standard_thermal_params[:k]
    @test body.point_params[1].kp == 6*standard_thermal_params[:k] / (π * standard_thermal_params[:δ]^4)
    @test k_linear(1.0, body.point_params[1].kp, 300.0) == body.point_params[1].kp
    @test k_linear(1.0, body.point_params[1].kp, 600.0) == 2*body.point_params[1].kp
    @test cv_linear(1.0, 300.0) == body.point_params[1].cv
    @test cv_linear(1.0, 600.0) == 2*body.point_params[1].cv
end