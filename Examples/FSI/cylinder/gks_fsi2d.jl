using HyperFSI

include("mesh.jl")

#--- flow setup---#
set = Setup(
    case = "rectangle",
    space = "2d0f0v",
    boundary = ["fix", "extra", "extra", "extra"],
    limiter = "vanleer",
    cfl = 0.2,
    maxTime = 2.0,
    flux = "gks",
    hasForce = false,
    )

#--- Dimensionless factor of flow variables ---#
ref = HyperFSI.ReferenceVariables(;
    length = 1.0,
    temperature = 294.0,
    velocity = 417.2881,
    density = 0.603,
    pressure = 52500 * 2,
    k=1.0
    )

nt = 1e5 # Number of time steps for coupling simulation
n_points = 50  
x_min, x_max = -2.0, -5.0 
t = range(0.0, 1.0+1.0/n_points, length=n_points+2)

a = 2.0  
X1 = reverse(x_min .+ (x_max - x_min) .* (exp.(a * t) .- 1.0) ./ (exp(a) - 1.0))

X2 = range(-2.0, 2.0, length=257)

n_points = 50 
x_min, x_max = 2.0, 5.0 

t = range(0.0, 1.0+1.0/n_points, length=n_points+2)

a = 2.0  
X3 = x_min .+ (x_max - x_min) .* (exp.(a * t) .- 1.0) ./ (exp(a) - 1.0)

X = unique(vcat(X1, X2, X3))

n_points = 50 
x_min, x_max = -2.0, -5.0 


t = range(0.0, 1.0+1.0/n_points, length=n_points+2)

a = 2.0 
Y1 = reverse(x_min .+ (x_max - x_min) .* (exp.(a * t) .- 1.0) ./ (exp(a) - 1.0))

Y2 = range(-2.0, 2.0, length=257)

n_points = 50 
x_min, x_max = 2.0, 5.0 

t = range(0.0, 1.0+1.0/n_points, length=n_points+2)
a = 2.0  
Y3 = x_min .+ (x_max - x_min) .* (exp.(a * t) .- 1.0) ./ (exp(a) - 1.0)

Y = unique(vcat(Y1, Y2, Y3))

ps = build_mesh(-5.0, 5.0, size(X, 1) - 3, -5.0, 5.0, size(Y, 1) - 3, X, Y)
vs = nothing 

U_∞ = 0.8 * sound_speed(1.0, 1.4)
Re = 2000.0

gas = Gas(Kn = 1e-3, Ma = 0.8, Pr = 0.71, K = 3.0, γ = 1.4, μᵣ = U_∞ / Re, ω = 0.74)

 # Initialization for flow variables; [density, velocity_x, velocity_y, pressure]
prim0 = [1.0, 0.0, 0.0, 1.0] # simulation domain
prim1 = [1.0, U_∞, 0.0, 1.0] # inflow condition

bc(x, y, args...) = x < 10.0 ? prim1 : prim0
fw(x, y, args...) = prim_conserve(bc(x, y, args...), gas.γ)

ib0 = IB(fw, bc, NamedTuple())

ks = SolverSet(set, ps, vs, gas, ib0)
ctr, a1face, a2face = init_fvm(ks; structarray = true)

dtf = timestep(ks, ctr, 0.0) * ref.time # stable timestep for flow simulation
fv = Flowstep(steps=nt, stepsize=dtf) #Time solver setup for flow simulation

# Structure setup
#Reading Mesh Geometry and Identifying Boundary Physical Groups
#Lines 5–8 define the outer boundary in the mesh.
geo = Post2D("2d_output.inp",["Line5","Line6","Line7","Line8"]) 

bc_index = geo.bc_nodes # Boundary node indices in PD node set
pos = geo.pos # Position of all the PD nodes
area = geo.area # Area of all the PD nodes
height = sqrt(minimum(area)) # Height of the cylinder for 2D
vol = height * area #Volume of all the PD nodes

body = Body(BBTMMaterial(), pos, vol)
material!(body; horizon=4height, E=2.0e9, rho=1e3, epsilon_c=1.0e-2, kc=11400.0, aph=1.0e-6, cv=1.0, rft=294.0, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=height) #thermomechanical parameters

bcst = Bcstruct(pos) #

hsource_databc!(body, bcst.hsource, :all_points) #Thermal flux boundary condition
Peridynamics.forcedensity_databc!(body, bcst.pressure, :all_points, [1,2,3])

vv = Thermomechstep(steps=nt) # Time solver setup for thermal diffusion simulation in structure

output_path = joinpath(@__DIR__, "cylinder")
rm(output_path;force=true,recursive=true)
mkdir(output_path)

# FSI job setup
job = FSI_job(ks, body, fv, vv; path=joinpath(@__DIR__, "cylinder"), freq=20, fields=(:displacement, :temperature, :damage))

# Submit the job to the PeriFlowX solver
# "T" = heat transfer mode (active),
# "TM" = thermomechanical interaction modes (active)
ret = FSI_submit(job, bcst, geo, "TM"; ref = ref)
