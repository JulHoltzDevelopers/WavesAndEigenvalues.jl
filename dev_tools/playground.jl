##
include("make_FEM_header.jl")
##
m, n, k = SymPy.symbols("m n k")
tet_integrate(x^m*y^n*z^k)
