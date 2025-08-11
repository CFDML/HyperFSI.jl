@testitem "Dualstep" begin
    using HyperFSI: Dualstep
    v0 = Dualstep(steps=1, ADRsteps=10) 
    v1 = Dualstep(steps=1, stepsize = 1e-6, ADRsteps=10) 
    v2 = Dualstep(time = 0.5, steps = 1, ADRsteps=10)
    v3 = Dualstep(time = 0.5, steps = 1, safety_factor = 0.5, ADRsteps=10)
    v4 = Dualstep(time = 0.5, steps = 1, safety_factor = 0.7, ADRsteps=100, d_factor=0.2)
    v5 = Dualstep(time = 0.5, steps = 1, safety_factor = 0.9, ADRsteps=1000, d_factor=0.2, ADRerror=1e-7)

    @test typeof(v0) == Dualstep
    @test v0.end_time == -1
    @test v0.n_steps == 1
    @test v0.Δt == -1
    @test v0.safety_factor == 0.7
    @test v0.ADRn_steps == 10
    @test v0.Λ == 1.0
    @test v0.ADRerror == 1e-5
    @test v1.end_time == -1
    @test v1.n_steps == 1
    @test v1.Δt == 1e-6
    @test v1.safety_factor == 0.7
    @test v1.ADRn_steps == 10
    @test v1.Λ == 1.0
    @test v1.ADRerror == 1e-5
    @test v2.end_time == 0.5
    @test v2.n_steps == 1
    @test v2.Δt == -1
    @test v2.safety_factor == 0.7
    @test v2.ADRn_steps == 10
    @test v2.Λ == 1.0
    @test v2.ADRerror == 1e-5
    @test v3.end_time == 0.5
    @test v3.n_steps == 1
    @test v3.Δt == -1
    @test v3.safety_factor == 0.5
    @test v3.ADRn_steps == 10
    @test v3.Λ == 1.0
    @test v3.ADRerror == 1e-5
    @test v4.end_time == 0.5
    @test v4.n_steps == 1
    @test v4.Δt == -1
    @test v4.safety_factor == 0.7       
    @test v4.ADRn_steps == 100
    @test v4.Λ == 0.2
    @test v4.ADRerror == 1e-5
    @test v5.end_time == 0.5    
    @test v5.n_steps == 1
    @test v5.Δt == -1
    @test v5.safety_factor == 0.9
    @test v5.ADRn_steps == 1000
    @test v5.Λ == 0.2
    @test v5.ADRerror == 1e-7
end