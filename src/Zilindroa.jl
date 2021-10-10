module Zilindroa

#Julia Libs
using LinearAlgebra
using LinearAlgebra.BLAS
using DifferentialEquations
using SparseArrays
using Dates
using Parameters
using HDF5
using DelimitedFiles
using ProgressMeter

#Import structures
include("./structs/Operators.jl")
include("./structs/Parameters.jl")

#Operators 
include("./operators/differential_operator.jl")
include("./operators/BC_operators.jl")

#Simulation helpers
include("./sim_utils/get_gaussian_pulse.jl")  

include("./sim_utils/setup_operators.jl")       
include("./sim_utils/setup_folders.jl")

include("./sim_utils/save_parameters_to_file.jl")

include("./sim_utils/print_parameters_to_screen.jl")

include("./sim_utils/time_to_seconds.jl")
include("./sim_utils/exit_message.jl")
include("./sim_utils/welcome_message.jl")

#Solvers
include("./solvers/Zilindroa_solver.jl")

end