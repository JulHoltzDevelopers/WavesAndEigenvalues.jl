import Pkg
Pkg.activate(".")
Pkg.instantiate()
cd("./examples/shape")

using WavesAndEigenvalues.Helmholtz
import WavesAndEigenvalues.Meshutils: insert_smplx!, find_smplx,sort_smplx
import Interpolations

include("read_vtu_wild.jl")
#mesh=octosplit(mesh)
##
#mesh=Mesh("ethv4.msh",scale=.001)
#collect_lines!(mesh)
#z=[x_min, 0.8, 0.80076, 0.811352, 0.819826, 0.856, 0.8685, 0.966, 0.978, 1.059, 1.0698, 1.16, 1.17, 1.18, 1.19, 1.27, 1.28, 1.34, 1.35, 1.44, 1.45, 1.52, 1.53, 1.67, 1.68, 1.87, 1.89, 2.00, 2.01, 2.04, 2.05, 2.07, 2.08, x_max]
#T=[347,  347, 714.2, 714.2, 676.5, 676.5, 813.6, 813.6, 839.5, 839.5, 745, 745, 740, 740, 718, 718, 719.8, 719.8, 700, 700, 702, 702, 811, 811, 838, 838, 839, 839, 790, 790, 837, 837, 839, 839]
#T = Interpolations.LinearInterpolation(z, T, extrapolation_bc = Interpolations.Flat())

z=[0 300;
 0.29837 300;
 0.33407 616.1928;
 0.357 1353.6472;
 0.37757596863661 1213.4379;
 0.457491304166493 1768.7575;
 0.517885848938566 1890.6961;
 0.555344663635985 1890.6961;
 0.587 1825.2656;
 0.624890478375109 1435.6567;
 0.642732577053009 1474.3202;
 0.701726613004129 1474.3202;
 0.741055970304876 1347.1389;
 0.770888726696417 1373.9059;
 0.795541518955666 1373.9059;
 0.807532176669308 1382.8283;
 0.842448971931434 1382.8283;
 0.85 1365]
T = Interpolations.LinearInterpolation(z[:,1], z[:,2], extrapolation_bc = Interpolations.Flat())
R=288.68 #gas cosntant of air
γ=1.4 #ration of specific heats
##
function sosfield(mesh,c_type=0)
    function speedofsound(x,y,z)
     #return T(x)
     return sqrt(γ*R*T(x-1.32+0.3775))
    end
    if c_type==0
        C=zeros(Float64,length(mesh.tetrahedra))#*347.0
        for idx=1:length(mesh.tetrahedra)
            cntr=sum(mesh.points[:,mesh.tetrahedra[idx]],dims=2)/4 #centerpoint of a tetrahedron is arithmetic mean of the vertices
            C[idx]=speedofsound(cntr...)
        end
    else
        C=zeros(Float64,size(mesh.points,2))
        for idx in 1:size(mesh.points,2)
            pnt=mesh.points[:,idx]
            C[idx]=speedofsound(pnt...)
        end
    end

    return C
end
C=sosfield(mesh,0)


##

dscrp=Dict()
dscrp["Interior"]=(:interior,())
#D["inlet"]=(:admittance, (:Y_in,1E15))
dscrp["Outlet"]=(:admittance, (:Y_out,1E15))
#D["damper"]=(:admittance, (:Y,1E15))
#D["Flame"]=(:flame,(γ,ρ,Q02U0,x_ref,n_ref,:n,:τ,1,.001))
D["Flame"]=(:fancyflame,(γ,ρ,Q02U0,x_ref,n_ref,[:n1,:n2],[:τ1,:τ2],[:a1,:a2],[1.0,.5],[.001,0.0005],[.01,.01]))

#unify!(mesh,"Flame2","Flame")
#D["Flame"]=(:fancyflame,(γ,ρ,Q02U0,x_ref,n_ref,:n1,:τ1,:a1,1.0,.001,.01,),)
#D["Flame2"]=(:fancyflame,(γ,ρ,Q02U0,x_ref,n_ref,:n1,:τ1,:a1,.5,.0005,.01,),)
L=discretize(mesh, dscrp, C, el_type=1, c_type=0)
#L.params[:n]=1#compute_volume!(mesh,"Flame")/5000
#L.params[:Y]=1E20

#param=:τ
##
sol,n,flag=householder(L,(150)*2*pi,maxiter=30,order=5,tol=1E-10,n_eig_val=1,output=true);

## shape_sensitivity
surface_points, tri_mask, tet_mask =get_surface_points(mesh)
normal_vectors=get_normal_vectors(mesh)
sens=discrete_adjoint_shape_sensitivity(mesh,dscrp,C,surface_points,tri_mask,tet_mask,L,sol)
normed_sens=normalize_sensitivity(surface_points,normal_vectors,tri_mask,sens)
normal_sens = normal_sensitivity(normal_vectors, normed_sens)

DD=Dict()
DD["normed_sens1"]=real.(normed_sens[1,:])
DD["normed_sens2"]=real.(normed_sens[2,:])
DD["normed_sens3"]=real.(normed_sens[3,:])
DD["normal_sens"]=real.(normal_sens)
vtk_write_tri("eth_tri",mesh,DD)

##

mode="[$(string(round((sol.params[:ω]/2/pi),digits=2)))]Hz"

data=Dict()
data["abs_p"*mode]=abs.(sol.v)./maximum(abs.(sol.v))
data["phase_p"*mode]=angle.(sol.v)
data["abs_p_adj"*mode]=abs.(sol.v_adj)./maximum(abs.(sol.v_adj))
data["phase_p_adj"*mode]=angle.(sol.v_adj)
data["speed_of_sound"]=C
data["discrete adjoint x"*mode]=real.(sens[1,:])
data["discrete adjoint y"*mode]=real.(sens[2,:])
data["discrete adjoint z"*mode]=real.(sens[3,:])
vtk_write("eth",mesh,data)
##
## FD approach
fd_sens=forward_finite_differences_shape_sensitivity(mesh,dscrp,C,surface_points,tri_mask,tet_mask,L,sol)
## write output
mode="[$(string(round((sol.params[:ω]/2/pi),digits=2)))]Hz"
data=Dict()
data["abs_p"*mode]=abs.(sol.v)./maximum(abs.(sol.v))
data["phase_p"*mode]=angle.(sol.v)
data["abs_p_adj"*mode]=abs.(sol.v_adj)./maximum(abs.(sol.v_adj))
data["phase_p_adj"*mode]=angle.(sol.v_adj)
data["discrete adjoint x"*mode]=real.(fd_sens[1,:])
data["discrete adjoint y"*mode]=real.(fd_sens[2,:])
data["discrete adjoint z"*mode]=real.(fd_sens[3,:])
vtk_write("eth_fd", mesh, data)

##
#(type Γ by typing \Gamma and Enter)
# Γ=[20.0+5.0im, 20.0-5.0im, 400.0-5.0im, 400.0+5.0im].*2*pi #corner points for the contour (in this case a rectangle)
# Ω, P = beyn(L,Γ,l=10,N=64, output=true)
# ##verify beyn's solution by local householder iterationtol=1E-10
# D=Dict()
# tol=1E-10
# for idx=1:length(Ω)
#     sol,nn,flag=householder(L,Ω[idx],maxiter=6,order=5,tol=tol,n_eig_val=3,output=false);
#     mode="[$(string(round((Ω[idx]/2/pi),digits=2)))]Hz"
#     if abs(sol.params[:ω]-Ω[idx])<1
#         println("Kept:$mode")
#         D["abs:"*mode]=abs.(sol.v)./maximum(abs.(sol.v)) #note that in Julia any function can be performed elementwise by putting a . before the paranthesis
#         D["phase:"*mode]=angle.(sol.v)
#     else
#         convmode="[$(string(round((sol.params[:ω]/2/pi),digits=2)))]Hz"
#         println("Excluded:$mode because converged to $convmode")
#     end
# end
# ##
# D["speedofsound"]=C
# vtk_write("ETH",mesh,D)
