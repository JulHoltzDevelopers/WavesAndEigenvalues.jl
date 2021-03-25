path=pwd()

##

using Documenter
import Pkg
Pkg.activate(".")
using WavesAndEigenvalues
##
cd("./docs")
makedocs(format = Documenter.HTML(
        prettyurls = false,
        ),
        sitename="WavesAndEigenvalues.jl",
        authors = "Georg A. Mensah, Alessandro Orchini, and contributors.",
        modules =[WavesAndEigenvalues],
        pages = [
        "Home" => "index.md",
        "Install" => "installation.md",
        "Examples" => Any["tutorial_00_NLEVP.md",
                        "tutorial_01_rijke_tube.md",
                        "tutorial_02_global_eigenvalue_solver.md",
                        "tutorial_03_local_eigenvalue_solver.md",
                        "tutorial_04_perturbation_theory.md",
                        "tutorial_05_mesh_refinement.md",
                        "tutorial_06_second_order_elements.md",
                        ],
        "API" => Any["NLEVP.md","Meshutils.md", "Helmholtz.md"]
        ]
        )
cd(path)

deploydocs(
    repo = "github.com/JulHoltzDevelopers/WavesAndEigenvalues.jl.git",
    push_preview=true
    )
