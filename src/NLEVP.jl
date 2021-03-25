"""
Module containing routines to solve and perturb nonlinear eigenvalue problems.
"""
module NLEVP
include("NLEVP_exports.jl")
#header
using SparseArrays

include("./NLEVP/algebra.jl")
include("./NLEVP/polys_pade.jl")
include("./NLEVP/LinOpFam.jl")
include("./NLEVP/perturbation.jl")
include("./NLEVP/Householder.jl")
include("./NLEVP/nicoud.jl")
include("./NLEVP/picard.jl")
include("./NLEVP/iterative_solvers.jl")
#include("./NLEVP/mehrmann.jl")
include("./NLEVP/beyn.jl")
include("./NLEVP/solver.jl")
include("./NLEVP/toml.jl")
include("./NLEVP/save.jl")
include("./NLEVP/gallery.jl")
#using .TOML
end
