using Documenter
using Pkg
Pkg.activate("./")
using WavesAndEigenvalues


makedocs(
    format = Documenter.HTML(
            prettyurls = false,
    ),
    sitename = "WavesAndEigenvalues",
    authors = "Georg A. Mensah, Alessandro Orchini, and contributors.",
    modules = [WavesAndEigenvalues],
    pages = [
    "Home" => "index.md",
    "Documentation" => Any[
        "Mesh.md",
        "Helmholtz.md",
        "NLEVP.md",
     ],
     "Examples" => Any[
     "load_mesh.md",
     "tutorial_01_rijke_tube.md",
     "tutorial_02_gloabal_eigenvalue_solver.md",
     "tutorial_03_local_eigenvalue_solver.md",
     ],
 ],
)

deploydocs(
    repo = "github.com/JulHoltzDevelopers/WavesAndEigenvalues.jl.git",
    push_preview=true
    )
