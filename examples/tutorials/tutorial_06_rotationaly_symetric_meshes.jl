## import Helmholtz-solver library
using WavesAndEigenvalues.Helmholtz
## load mesh from file
mesh=Mesh("NTNU.msh",scale=1.0)
#create unit and full mesh from half-cell
doms=[("Interior",:full),("Inlet",:full), ("Outlet_high",:full), ("Outlet_low",:full), ("Flame",:unit),]#
collect_lines!(mesh)
unit_mesh=extend_mesh(mesh,doms,unit=true)
full_mesh=extend_mesh(mesh,doms,unit=false)

## Helmholtz solver
#define speed of sound field
function speedofsound(x,y,z)
    if z<0.415
        return 347.0#m/s
    else
        return 850.0#m/s
    end
end


c=generate_field(unit_mesh,speedofsound,0)
C=generate_field(full_mesh,speedofsound,0)
#D=Dict()

#vtk_write("NTNU_c",full_mesh,D)
##describe model
D=Dict()
D["Interior"]=(:interior,())
D["Outlet_high"]=(:admittance, (:Y_in,1E15))
D["Outlet_low"]=(:admittance, (:Y_out,1E15))
##discretize models
L=discretize(full_mesh, D, C, el_type=1, c_type=0)
l=discretize(unit_mesh, D, c, el_type=1, c_type=0,b=:b)
## solve model using beyn
#(type Γ by typing \Gamma and Enter)
Γ=[150.0+5.0im, 150.0-5.0im, 1000.0-5.0im, 1000.0+5.0im].*2*pi #corner points for the contour (in this case a rectangle)
Ω, P = beyn(L,Γ,l=10,N=64, output=true)
ω, p = beyn(l,Γ,l=10,N=128, output=true)
##verify beyn's solution by local householder iteration
tol=1E-10
D=Dict()
for idx=1:length(Ω)
    sol,nn,flag=householder(L,Ω[idx],maxiter=6,order=5,tol=tol,n_eig_val=3,output=false);
    mode="[$(string(round((Ω[idx]/2/pi),digits=2)))]Hz"
    if abs(sol.params[:ω]-Ω[idx])<5*tol
        println("Kept:$mode")
        D["abs:"*mode]=abs.(sol.v)./maximum(abs.(sol.v)) #note that in Julia any function can be performed elementwise by putting a . before the paranthesis
        D["phase:"*mode]=angle.(sol.v)
    else
        convmode="[$(string(round((sol.params[:ω]/2/pi),digits=2)))]Hz"
        println("Excluded:$mode because converged to $convmode")
    end
end
##
D=Dict()
D["speedofsound"]=C
vtk_write("NTNU",full_mesh,D)
##

sol,nn,flag=householder(l,1000*2*pi,maxiter=6,order=1,tol=1E-9,n_eig_val=1,output=true);
mode="[$(string(round((sol.params[:ω]/2/pi),digits=2)))]Hz"
v=bloch_expand(unit_mesh,sol)
D["abs Bloch:"*mode]=abs.(v)./maximum(abs.(v)) #note that in Julia any function can be performed elementwise by putting a . before the paranthesis
D["phase Bloch:"*mode]=angle.(v)
#ol,nn,flag=householder(L,sol.params[:ω],maxiter=6,order=5,tol=1E-9,n_eig_val=3,output=true);
Sol,nn,flag=householder(L,sol.params[:ω],maxiter=6,order=5,tol=1E-5,n_eig_val=3,output=true);

mode="[$(string(round((Sol.params[:ω]/2/pi),digits=2)))]Hz"
D["abs:"*mode]=abs.(Sol.v)./maximum(abs.(Sol.v)) #note that in Julia any function can be performed elementwise by putting a . before the paranthesis
D["phase:"*mode]=angle.(Sol.v)
vtk_write("NTNU",full_mesh,D)

##

##

D=Dict()
D["speedofsound"]=C
Sol,nn,flag=householder(L,500*2*pi,maxiter=6,order=5,tol=1E-5,n_eig_val=3,output=true);
mode="[$(string(round((Sol.params[:ω]/2/pi),digits=2)))]Hz"
D["abs:"*mode]=abs.(Sol.v)./maximum(abs.(Sol.v)) #note that in Julia any function can be performed elementwise by putting a . before the paranthesis
D["phase:"*mode]=angle.(Sol.v)
vtk_write("NTNU",full_mesh,D)

##
