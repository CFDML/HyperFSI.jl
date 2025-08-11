@testitem "Thermomechstep" begin
    using HyperFSI: Thermomechstep
    v0 = Thermomechstep(steps=1) 
    v1 = Thermomechstep(steps=1, stepsize = 1e-6) 
    v2 = Thermomechstep(time = 0.5, steps = 1)
    v3 = Thermomechstep(time = 0.5, steps = 1, safety_factor = 0.5)

    @test typeof(v0) == Thermomechstep
    @test v0.end_time == -1
    @test v0.n_steps == 1
    @test v0.Δt == -1
    @test v0.safety_factor == 0.7
    @test v0.Δt_mech == -1
    @test v0.Δt_therm == -1
    @test v1.end_time == -1
    @test v1.n_steps == 1
    @test v1.Δt == 1e-6
    @test v1.safety_factor == 0.7
    @test v1.Δt_mech == -1
    @test v1.Δt_therm == -1
    @test v2.end_time == 0.5
    @test v2.n_steps == 1
    @test v2.Δt == -1
    @test v2.safety_factor == 0.7
    @test v2.Δt_mech == -1
    @test v2.Δt_therm == -1
    @test v3.end_time == 0.5
    @test v3.n_steps == 1
    @test v3.Δt == -1
    @test v3.safety_factor == 0.5
    @test v3.Δt_mech == -1
    @test v3.Δt_therm == -1
end