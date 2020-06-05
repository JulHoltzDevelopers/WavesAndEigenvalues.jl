using Documenter
#AO: I need these two lines to build docs for the LOCAL version
using Pkg
Pkg.activate("../../WavesAndEigenvalues_devel/")
using WavesAndEigenvalues

push!(LOAD_PATH,"../src/")

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
     ],
 ],
)
