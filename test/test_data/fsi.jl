using HyperFSI

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
    geo = Post2D("Ring.inp", ["Line5","Line6","Line7","Line8"])  
    bc_index = geo.bc_nodes  # Boundary node indices
    pos = geo.pos            # Node positions [m]
    area = geo.area          # Nodal areas [m²]
    height = sqrt(minimum(area))  # Characteristic length
    vol = height * area      # Nodal volumes [m³]
    body = Peridynamics.Body(BBTMaterial(), pos, vol)

    Peridynamics.material!(body;
        horizon = 4*height,    # Peridynamic horizon (4× characteristic length)
        E = 2.0e11,            # Elastic modulus
        rho = 1e3,             # Material density
        epsilon_c = 1e-2,      # Critical fracture strain
        kc = 11400.0,          # Thermal conductivity
        aph = 1.0e-6,          # Thermal expansion
        cv = 1.0,              # Heat capacity
        rft = 273,             # Reference temperature [K]
        h = 25.0,              # Convection coefficient [W/(m²·K)]
        hσ = 5.73e-8,          # Stefan-Boltzmann Constant 
        hϵ = 0.9,              # Emissivity
        tem∞ = 0.0,            # Ambient temperature [K]
        thick = height         # Thickness for 2D
    )
    bcst = Bcstruct(pos)  

    hsource_databc!(body, bcst.hsource, :all_points)  
    vv = Thermstep(steps=1)  
    job = FSI_job(ks, body, fv, vv;
        path = joinpath(@__DIR__, "cylinder"), # - Output directory "./cylinder"
        freq = 1,                              # Output frequency
        fields = (:displacement, :temperature, :damage)# - Output fields
    )
    ret = FSI_submit(job, bcst, geo, "T"; ref = ref)
