@testitem "FSI_job" begin
    using HyperFSI: FSI_job
        
    set = KitBase.Setup(
        case = "rectangle",       # Case identifier
        space = "2d0f0v",         # 2D simulation space
        boundary = ["fix", "extra", "extra", "extra"],  # Boundary conditions
        limiter = "vanleer",       # Flux limiter type
        cfl = 0.2,                # CFL number for stability
        maxTime = 2.0,            # Maximum simulation time
        flux = "gks",             # Flux calculation scheme
        hasForce = false,         # External forces flag
    )
    
    ref = HyperFSI.ReferenceVariables(;
        length = 1.0,            # Reference length [m]
        temperature = 294.0,      # Reference temperature [K]
        velocity = 417.2881,      # Reference velocity [m/s] !! Must be Sqrt(2RT)!!
        density = 0.603,          # Reference density [kg/m³]
        pressure = 52500 * 2,     # Reference pressure [Pa]
        k = 1.0                   # Thermal diffusion factor
    )

    ps = KitBase.PSpace2D(-3.0, 3.0, 140, -2.0, 2.0, 90, 1, 1)
    vs = nothing

    U_∞ = 0.8 * KitBase.sound_speed(1.0, 1.4)  # Freestream velocity
    Re = 2000.0                                # Reynolds number
    gas = KitBase.Gas(
        Kn = 1e-3,       # Knudsen number
        Ma = 0.8,        # Mach number  
        Pr = 0.71,       # Prandtl number
        K = 3.0,         # Internal DOF!!!
        γ = 1.4,         # Specific heat ratio
        μᵣ = U_∞ / Re,   # Reference viscosity
        ω = 0.74         # Viscosity exponent
    )
    prim0 = [1.0, 0.0, 0.0, 1.0]
    prim1 = [1.0, U_∞, 0.0, 1.0] 
    bc(x, y, args...) = x < -1.0 ? prim1 : prim0
    fw(x, y, args...) = KitBase.prim_conserve(bc(x, y, args...), gas.γ)

    ib0 = KitBase.IB(fw, bc, NamedTuple())
    ks = KitBase.SolverSet(set, ps, vs, gas, ib0)
    ctr, a1face, a2face = KitBase.init_fvm(ks; structarray = true)
    dtf = KitBase.timestep(ks, ctr, 0.0) * ref.time  
    fv = Flowstep(steps=1, stepsize=dtf) 

    mesh_file = "/Users/shiweihu/Codes/Develops/HyperFSI/test/test_data/Ring.inp"
    geo = Post2D(mesh_file, ["Line5","Line6","Line7","Line8"])  
    pos = geo.pos            # Node positions [m]
    area = geo.area          # Nodal areas [m²]
    height = sqrt(minimum(area))  # Characteristic length
    vol = height * area      # Nodal volumes [m³]
    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=4height, E=2.0e11, rho=1e3, epsilon_c=1e-2, kc=11400.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height)
    bcst = Bcstruct(pos) #

    hsource_databc!(body, bcst.hsource, :all_points)
    vv = Thermstep(steps=1)  
    job = FSI_job(ks, body, fv, vv;
        path = joinpath(@__DIR__, "cylinder"), # - Output directory "./cylinder"
        freq = 1,                              # Output frequency
        fields = (:displacement, :temperature, :damage)# - Output fields
    )

    @test job isa FSI_job
    @test job.f_time_solver isa Flowstep
    @test job.s_time_solver isa Thermstep
    @test job.spatial_setup isa Peridynamics.Body
    @test job.flow_setup isa KitBase.SolverSet
    @test job.options.export_allowed == true
    @test job.options.freq == 1
    @test job.options.root == joinpath(@__DIR__, "cylinder")
    @test job.options.vtk == joinpath(@__DIR__, "cylinder/structure")
    @test job.options.fields == [:displacement, :temperature, :damage]
    @test job.options.flow == joinpath(@__DIR__, "cylinder/flow")
end

@testitem "Bcstruct" begin
    using HyperFSI: Bcstruct

    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]

    bc = Bcstruct(pos) 
    @test bc isa Bcstruct
    @test size(bc.hsource) == (1, 4)
    @test size(bc.pressure) == (3, 4)
    @test all(bc.hsource .== 0.0)
    @test all(bc.pressure .== 0.0)
end