"""
Module containing routines to solve and perturb nonlinear eigenvalue problems.
"""
module NLEVP
include("NLEVP_exports.jl")
#header
using SparseArrays

include("./algebra.jl")
include("./LinOpFam.jl")
include("./perturbation.jl")
include("./Householder.jl")
include("./beyn.jl")
include("./toml.jl")
include("./save.jl")
#using .TOML
end
