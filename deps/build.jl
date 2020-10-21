include("../src/perturbation.jl")
println("Building WavesAndEigenvalues.jl...")

if "JULIA_WAE_PERT_ORDER" in keys(ENV)
    N=ENV["JULIA_WAE_PERT_ORDER"]
    N=parse(Int,N)
else
    N=10
end
println(" for order $N...")
disk_MN(N)
println("Done!")
