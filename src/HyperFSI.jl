module HyperFSI

using MPI
using Reexport
@reexport using Peridynamics
@reexport using KitBase
using Base.Threads
using TimerOutputs
using ProgressMeter
using PointNeighbors
using LinearAlgebra
using DelimitedFiles
using Glob
using WriteVTK
using Dates

abstract type AbstractFlowTimeSolver end
abstract type AbstractPDGeometry end

# Pre processing
export Post2D, Post3D

# Coupling with DSMC_sparta
#export glob, write_sprata_bc_file_2d, read_bc_from_sparta, write_sparta_files, writedlm, readdlm

# New Material models
export BBTMaterial, BBTMMaterial, BBTTMaterial, BBTTAMaterial, BBTMAMaterial, BBTxyMaterial, BBTTxyMaterial


# New Discretization
export hsource_bc!, hsource_databc!, temperature_ic!, temperature_bc!, temperature_databc!, second_bcs!, find_sec_bcs_points 

# Running simulations
export Thermstep, Thermomechstep, Dualstep, Thermstep_ablation, Thermomechstep_ablation, Dualstep_ablation, Thermstep_anis, 
        Flowstep, FSI_job, FSI_submit, IBM2D, Bcstruct


include("IBM/boundary_counter.jl")
include("IBM/post_2d.jl")
include("IBM/post_3d.jl")
include("IBM/ibm.jl")
include("IBM/find_new_bc_edges.jl")

include("structure/physics/modify.jl")
include("structure/physics/new_boundaries.jl")
include("structure/physics/ablation.jl")

include("structure/time_solvers/thermstep.jl")
include("structure/time_solvers/thermomechstep.jl")
include("structure/time_solvers/dual_timesteps.jl")
include("structure/time_solvers/thermstep_ablation.jl")
include("structure/time_solvers/thermomechstep_ablation.jl")
include("structure/time_solvers/dual_timesteps_ablation.jl")
include("structure/time_solvers/thermstep_anis.jl")

include("structure/physics/bond_based_thermal_diffusion.jl")
include("structure/physics/bond_based_thermomechanics.jl")
include("structure/physics/bond_based_dualstep_thermomechanics.jl")
include("structure/physics/bond_based_thermal_diffusion_temp_dependent.jl")
include("structure/physics/bond_based_thermal_diffusion_ablation.jl")
include("structure/physics/bond_based_thermomechanics_ablation.jl")
include("structure/physics/bond_based_dualstep_thermomechanics_ablation.jl")

include("structure/physics/anisotropic/bond_based_thermal_diffusion.jl")
include("structure/physics/anisotropic/bond_based_thermal_diffusion_temp_dependent.jl")

include("fluid/FlowTimesolver.jl")
include("fluid/Evolution.jl")
include("fluid/Advance.jl")
include("fluid/output_vtk.jl")

#coupling_solvers
include("coupling_solvers/logs.jl") 
include("coupling_solvers/FSI_job.jl")
include("coupling_solvers/FSI_submit.jl") # central of jobs
end




