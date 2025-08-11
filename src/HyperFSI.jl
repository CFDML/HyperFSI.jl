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

# Pre processing
export Post2D, Post3D

# Coupling with DSMC_sparta
#export glob, write_sprata_bc_file_2d, read_bc_from_sparta, write_sparta_files, writedlm, readdlm

# New Material models
export BBTMaterial, BBTMMaterial


# New Discretization
export hsource_bc!, hsource_databc!, temperature_ic!, temperature_bc!, temperature_databc!, second_bcs!, find_sec_bcs_points 

# Running simulations
export Thermstep, Thermomechstep, Dualstep, Flowstep, FSI_job, FSI_submit, IBM2D, Bcstruct


include("IBM/post_2d.jl")
include("IBM/post_3d.jl")
include("IBM/ibm.jl")

include("structure/physics/modify.jl")
include("structure/physics/new_boundaries.jl")
include("structure/time_solvers/thermstep.jl")
include("structure/time_solvers/thermomechstep.jl")
include("structure/time_solvers/dual_timesteps.jl")
include("structure/physics/bond_based_thermal_diffusion.jl")
include("structure/physics/bond_based_thermomechanics.jl")
include("structure/physics/bond_based_dualstep_thermomechanics.jl")

include("fluid/FlowTimesolver.jl")
include("fluid/Evolution.jl")
include("fluid/Advance.jl")
include("fluid/output_vtk.jl")

#coupling_solvers
include("coupling_solvers/logs.jl") 
include("coupling_solvers/FSI_job.jl")
include("coupling_solvers/FSI_submit.jl") # central of jobs
end




