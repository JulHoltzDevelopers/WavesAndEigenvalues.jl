#import WavesAndEigenvalues.Meshutils: create_rotation_matrix_around_axis, get_rotated_index, get_reflected_index
using WavesAndEigenvalues.Helmholtz
## load mesh from file
case="NTNU"
mesh=Mesh("NTNU_12.msh",scale=1.0)
#create unit and full mesh from half-cell
doms=[("Interior",:full),("Inlet",:full), ("Outlet_high",:full), ("Outlet_low",:full), ("Walls",:full),("Flame",:unit),]#
collect_lines!(mesh)
unit_mesh=extend_mesh(mesh,doms,unit=true)
full_mesh=extend_mesh(mesh,doms,unit=false)
## meshing
#naxis=unit_mesh.dos.naxis
#nxbloch=unit_mesh.dos.nxbloch

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
dscrp=Dict()
dscrp["Interior"]=(:interior,())
dscrp["Outlet_high"]=(:admittance, (:Y_in,1E15))
dscrp["Outlet_low"]=(:admittance, (:Y_out,1E15))
##discretize models
L=discretize(full_mesh, dscrp, C, order=:1)
l=discretize(unit_mesh, dscrp, c, order=:1,b=:b)
## solve model using beyn
#(type Γ by typing \Gamma and Enter)
Γ=[150.0+5.0im, 150.0-5.0im, 1000.0-5.0im, 1000.0+5.0im].*2*pi #corner points for the contour (in this case a rectangle)
#Ω, P = beyn(L,Γ,l=10,N=64, output=true)
l.params[:b]=1
#ω, p = beyn(l,Γ,l=10,N=4*128, output=true)
##
sol, n, flag = householder(l,914*2*pi,maxiter=16, tol=0*1E-10,output = true, n_eig_val=3)




##
# identify surface points
surface_points, tri_mask, tet_mask =get_surface_points(unit_mesh)
# compute normal vectors on surface triangles
normal_vectors=get_normal_vectors(unit_mesh)
## DA approach
blochify_surface_points!(unit_mesh, surface_points, tri_mask, tet_mask)
sens=discrete_adjoint_shape_sensitivity(unit_mesh,dscrp,c,surface_points,tri_mask,tet_mask,l,sol)
##
# normalize point sensitivity with directed surface triangle area
normed_sens=normalize_sensitivity(surface_points,normal_vectors,tri_mask,sens)
# boundary-mass-based-normalizations
nsens = bound_mass_normalize(surface_points,normal_vectors,tri_mask,unit_mesh,sens)
# compute sensitivity in unit normal direction
normal_sens = normal_sensitivity(normal_vectors, normed_sens)

## correct sens for bloch formalism
#sens[:,unit_mesh.dos.naxis+1:unit_mesh.dos.naxis+unit_mesh.dos.nxbloch].+=sens[:,end-unit_mesh.dos.nxbloch+1:end]
sens[:,end-unit_mesh.dos.nxbloch+1:end]=sens[:,unit_mesh.dos.naxis+1:unit_mesh.dos.naxis+unit_mesh.dos.nxbloch]
## save data
DD=Dict()
DD["real normed_sens1"]=real.(normed_sens[1,:])
DD["real normed_sens2"]=real.(normed_sens[2,:])
DD["real normed_sens3"]=real.(normed_sens[3,:])
DD["real normal_sens"]=real.(normal_sens)
DD["imag normed_sens1"]=imag.(normed_sens[1,:])
DD["imag normed_sens2"]=imag.(normed_sens[2,:])
DD["imag normed_sens3"]=imag.(normed_sens[3,:])
DD["imag normal_sens"]=imag.(normal_sens)
vtk_write_tri(case*"_tri",unit_mesh,DD)


##
mode="[$(string(round((sol.params[:ω]/2/pi),digits=2)))]Hz"
data=Dict()
data["real DA x"*mode]=real.(sens[1,:])
data["real DA y"*mode]=real.(sens[2,:])
data["real DA z"*mode]=real.(sens[3,:])
data["imag DA x"*mode]=imag.(sens[1,:])
data["imag DA y"*mode]=imag.(sens[2,:])
data["imag DA z"*mode]=imag.(sens[3,:])
data["real nDA x"*mode]=real.(nsens[1,:])
data["real nDA y"*mode]=real.(nsens[2,:])
data["real nDA z"*mode]=real.(nsens[3,:])
data["imag nDA x"*mode]=imag.(nsens[1,:])
data["imag nDA y"*mode]=imag.(nsens[2,:])
data["imag nDA z"*mode]=imag.(nsens[3,:])
#data["finite difference"*mode]=real.(sens_FD)
#data["central finite difference"]=real.(sens_CFD)
#data["diff DA -- FD"*mode]=abs.((sens.-sens_FD))
vtk_write(case, unit_mesh, data)
data=Dict()
v=bloch_expand(unit_mesh,sol,:b)
data["abs_p"*mode]=abs.(v)./maximum(abs.(v))
data["phase_p"*mode]=angle.(v)
#data["abs_p_adj"*mode]=abs.(sol.v_adj)./maximum(abs.(sol.v_adj))
#data["phase_p_adj"*mode]=angle.(sol.v_adj)
vtk_write(case*"_modes", full_mesh, data)

## FD approach

fd_sens=forward_finite_differences_shape_sensitivity(unit_mesh,dscrp,c,surface_points,tri_mask,tet_mask,l,sol)
fd_sens[:,end-unit_mesh.dos.nxbloch+1:end]=fd_sens[:,unit_mesh.dos.naxis+1:unit_mesh.dos.naxis+unit_mesh.dos.nxbloch]


## write output
[]mode="[$(string(round((sol.params[:ω]/2/pi),digits=2)))]Hz"
data=Dict()
data["abs_p"*mode]=abs.(sol.v)./maximum(abs.(sol.v))
data["phase_p"*mode]=angle.(sol.v)
data["abs_p_adj"*mode]=abs.(sol.v_adj)./maximum(abs.(sol.v_adj))
data["phase_p_adj"*mode]=angle.(sol.v_adj)
data["real DA x"*mode]=real.(fd_sens[1,:])
data["real DA y"*mode]=real.(fd_sens[2,:])
data["real DA z"*mode]=real.(fd_sens[3,:])
data["imag DA x"*mode]=imag.(fd_sens[1,:])
data["imag DA y"*mode]=imag.(fd_sens[2,:])
data["imag DA z"*mode]=imag.(fd_sens[3,:])
#data["finite difference"*mode]=real.(fd_sens_FD)
#data["central finite difference"]=real.(sens_CFD)
#data["diff DA -- FD"*mode]=abs.((sens.-sens_FD))
vtk_write(case*"_fd", unit_mesh, data)
##
findmax(abs.(sens-fd_sens))
