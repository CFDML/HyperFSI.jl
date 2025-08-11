# HyperFSI.jl


## Overview
This computational framework is specifically designed for high-fidelity modeling of strongly coupled interactions between hypersonic flows and deformable structures with fracture mechanics. Its core capabilities include:

- **High-speed fluid-structure coupling**: Resolution of strongly coupled interactions between compressible hypersonic flows and nonlinear structural responses
- **Thermomechanical failure analysis**: Prediction of structural damage evolution (including complex fracture patterns) under extreme aerodynamic heating and dynamic loading conditions
- **Two-way multiphysics interactions**: Quantification of flow field modifications caused by structural topology changes

## Foundation
- **Flow solver**: [KitBase.jl](https://github.com/vavrines/KitBase.jl) (By Tianbai Xiao)
- **Structural solver**: Extended from [Peridynamics.jl](https://github.com/kaipartmann/Peridynamics.jl) (By Kai Partmann)
- **Coupling scheme**: Ghost-cell immersed boundary method (GCIBM)

## Key Capabilities
1. **Tightly coupled multiphysics solvers**:
   - Fluid-structure-thermal conduction analysis
   - Fluid-structure-mechanical analysis
   - Fully coupled **fluid-structure-thermomechanical** analysis
2. **Advanced structural modeling**:
   - Robust complex geometry handling
   - Thermomechanical failure prediction

## Ongoing Development
- Extension to **3D complex geometries**
- Extension to **more coupled fields** 
- GPU acceleration support

---

## Installation
```julia
julia> ] add HyperFSI

or

julia> using Pkg; Pkg.add("HyperFSI")

using HyperFSI
```


## Example usage
### Subsonic Flow Around a Deformable Cylinder with Thermomechanical Coupling

This example demonstrates a fully coupled fluid-structure-thermal simulation of subsonic flow interacting with a deformable cylinder structure. The simulation includes:
- Compressible flow solving using GKS (Gas-Kinetic Scheme)
- Peridynamics-based structural mechanics
- Two-way thermal-mechanical coupling
- Fracture prediction under aerodynamic heating

```julia
using HyperFSI

#--------------------------------------------------
# Flow Solver Configuration
#--------------------------------------------------
# Case setup for 2D rectangular domain with:
# - Left boundary: fixed inflow condition
# - Other boundaries: extra (non-reflective)
# - GKS flux scheme for compressible flow
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

#--------------------------------------------------
# Reference Variables (Non-dimensionalization)
#--------------------------------------------------
ref = HyperFSI.ReferenceVariables(;
    length = 1.0,            # Reference length [m]
    temperature = 294.0,      # Reference temperature [K]
    velocity = 417.2881,      # Reference velocity [m/s] !! Must be Sqrt(2RT)!!
    density = 0.603,          # Reference density [kg/m³]
    pressure = 52500 * 2,     # Reference pressure [Pa]
    k = 1.0                   # Thermal diffusion factor
)

#--------------------------------------------------
# Flow Domain Initialization
#--------------------------------------------------
# Physical space discretization:
# - x ∈ [-3,3] with 140 cells
# - y ∈ [-2,2] with 90 cells
ps = KitBase.PSpace2D(-3.0, 3.0, 140, -2.0, 2.0, 90, 1, 1)
vs = nothing

# Gas properties for air:
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

# Initial conditions [density, u, v, pressure]
# prim0: Domain initial condition 
# prim1: Inflow condition (Mach 0.8)
prim0 = [1.0, 0.0, 0.0, 1.0]
prim1 = [1.0, U_∞, 0.0, 1.0] 

# Boundary condition function:
# - Left boundary (x < -1.0): inflow condition
# - Elsewhere: domain initial condition
bc(x, y, args...) = x < -1.0 ? prim1 : prim0
fw(x, y, args...) = KitBase.prim_conserve(bc(x, y, args...), gas.γ)

# Initialize immersed boundary condition
ib0 = KitBase.IB(fw, bc, NamedTuple())

# Complete flow solver initialization
ks = KitBase.SolverSet(set, ps, vs, gas, ib0)
ctr, a1face, a2face = KitBase.init_fvm(ks; structarray = true)

# Time step calculation based on CFL condition
dtf = KitBase.timestep(ks, ctr, 0.0) * ref.time  
fv = Flowstep(steps=1, stepsize=dtf)  # Flow time integration setup

#--------------------------------------------------
# Structure Configuration (Peridynamics)
#--------------------------------------------------
# Load cylinder geometry from INP file:
# - "Ring.inp" contains cylinder mesh
# - Lines 5-8 define the outer boundary 
geo = Post2D("Ring.inp", ["Line5","Line6","Line7","Line8"])  

# Extract structural properties:
bc_index = geo.bc_nodes  # Boundary node indices
pos = geo.pos            # Node positions [m]
area = geo.area          # Nodal areas [m²]
height = sqrt(minimum(area))  # Characteristic length
vol = height * area      # Nodal volumes [m³]

# Create peridynamic body with Bond-based-thermal material model:
# - Bond-Based Thermomechanical material
# - Includes thermal expansion and heat conduction
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

# Initialize boundary condition for structure
bcst = Bcstruct(pos)  

# Apply thermal flux boundary condition 
hsource_databc!(body, bcst.hsource, :all_points)  

# Thermal diffusion solver setup
vv = Thermstep(steps=1)  

#--------------------------------------------------
# Coupled Simulation Setup
#--------------------------------------------------
# Create FSI job with:
job = FSI_job(ks, body, fv, vv;
    path = joinpath(@__DIR__, "cylinder"), # - Output directory "./cylinder"
    freq = 1,                              # Output frequency
    fields = (:displacement, :temperature, :damage)# - Output fields
)

# Submit simulation with "T" mode（Active heat transfer coupling）
ret = FSI_submit(job, bcst, geo, "T"; ref = ref)
```

## Contributing
We warmly welcome contributions of all kinds! Here's how you can contribute:

- Fork the repository
- Create a new branch (git checkout -b my-feature)
- Commit your changes (git commit -am 'Add some feature')
- Push to the branch (git push origin my-feature)
- Submit a pull request
- Contact: hushiwei@imech.ac.cn

## Authors
* [Shiwei Hu](https://www.researchgate.net/profile/Shiwei-Hu-7); hushiwei@imech.ac.cn
* [Mingshuo Han](); hanmingshuo@imech.ac.cn
* [Tianbai Xiao](https://people.ucas.edu.cn/~txiao); txiao@imech.ac.cn
* [Yonghao Zhang](https://people.ucas.ac.cn/~0075739); yonghao.zhang@imech.ac.cn

## Affiliation
[Centre for Interdisciplinary Research in Fluids, Institute of Mechanics, Chinese Academy of Sciences](https://imech.cas.cn/cirf/zxjj/)

## Acknowledgments:
We gratefully acknowledge Kai Partmann for his sustained technical support and insightful discussions.

## License
HyperFSI.jl is licensed under the MIT License.

## Citation
@misc{fsi,
  title={A Navier-Stokes-Peridynamics hybrid algorithm for the coupling of compressible flows and fracturing materials}, 
  author={Mingshuo Han and Shiwei Hu and Tianbai Xiao and Yonghao Zhang},
  year={2025},
  eprint={2504.11006},
  archivePrefix={arXiv},
  primaryClass={physics.comp-ph},
  url={[https://arxiv.org/abs/2504.11006](https://arxiv.org/abs/2504.11006)}
}

