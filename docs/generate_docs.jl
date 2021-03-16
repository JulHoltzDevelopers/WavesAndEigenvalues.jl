using Revise #src
Pkg.activate(".") #src
##
import Literate
##
path=pwd()
outputdir=path*"/docs/src"

exmpls=[("00_NLEVP",false),
        ("01_rijke_tube",false),
        ("02_global_eigenvalue_solver",false),
        ("03_local_eigenvalue_solver",false),
        ("04_perturbation_theory",true),
        ("05_mesh_refinement",false),
        ("06_second_order_elements",false)]

cd("./examples/tutorials")
for (exmp,exec) in exmpls
    Literate.markdown("./tutorial_$exmp.jl",outputdir,execute=true)
    #Literate.notebook("./tutorial_$exmp.jl",outputdir)
    #Literate.script("./tutorial_$exmp.jl",outputdir,keep_comments=true)
end

cd(path)
