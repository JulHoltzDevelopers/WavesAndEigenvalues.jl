#Simple tube example.
using WavesAndEigenvalues.Helmholtz
#import WavesAndEigenvalues.Meshutils: insert_smplx!, find_smplx
#import SparseArrays
#import ProgressMeter
#import LinearAlgebra
##
mesh=Mesh("Rijke_mm.msh",scale=0.001)
#mesh=Mesh("rijke_fine.msh",scale=0.01*0.5)
C=ones(length(mesh.tetrahedra))*347.0
dscrp=Dict()
dscrp["Interior"]=(:interior,())
dscrp["Outlet"]=(:admittance, (:Y,1E15))
##
L=discretize(mesh, dscrp, C)

##
sol, n, flag = householder(L,3*147*2*pi,maxiter=14, tol=1E-10,output = true, n_eig_val=3)




##
surface_points, tri_mask, tet_mask =get_surface_points(mesh)
normal_vectors=get_normal_vectors(mesh)

## DA approach


##
sens=discrete_adjoint_shape_sensitivity(mesh,dscrp,C,surface_points,tri_mask,tet_mask,L,sol)


##
normed_sens=normalize_sensitivity(surface_points,normal_vectors,tri_mask,sens)
##
normal_sens = normal_sensitivity(normal_vectors, normed_sens)
ω0=sol.params[:ω]
mode="[$(string(round((ω0/2/pi),digits=2)))]Hz"


##
DD=Dict()
DD["normed_sens1"]=real.(normed_sens[1,:])
DD["normed_sens2"]=real.(normed_sens[2,:])
DD["normed_sens3"]=real.(normed_sens[3,:])
DD["normal_sens"]=real.(normal_sens)
vtk_write_tri("tri",mesh,DD)


##
mode="[$(string(round((sol.params[:ω]/2/pi),digits=2)))]Hz"
data=Dict()
data["abs_p"*mode]=abs.(sol.v)./maximum(abs.(sol.v))
data["phase_p"*mode]=angle.(sol.v)
data["abs_p_adj"*mode]=abs.(sol.v_adj)./maximum(abs.(sol.v_adj))
data["phase_p_adj"*mode]=angle.(sol.v_adj)
data["discrete adjoint x"*mode]=real.(sens[1,:])
data["discrete adjoint y"*mode]=real.(sens[2,:])
data["discrete adjoint z"*mode]=real.(sens[3,:])
#data["finite difference"*mode]=real.(sens_FD)
#data["central finite difference"]=real.(sens_CFD)
#data["diff DA -- FD"*mode]=abs.((sens.-sens_FD))
vtk_write("shape_sense_v1", mesh, data)
## FD approach
sens_FD=forward_finite_differences_shape_sensitivity(mesh::Mesh,dscrp,C,surface_points,tri_mask,tet_mask,L,sol,h=1E-10)

## old FD implementation
# sens_FD=zeros(ComplexF64,size(mesh.points,2))
# p = ProgressMeter.Progress(length(surface_points),desc="FD Sensitivity... ", dt=1,
#      barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
#      barlen=20)
# for idx=1:length(surface_points)
#     #global G
#     pnt_idx=surface_points[idx]
#     crdnt=3
#     #forward finite difference
#     mesh_h=Mesh("mesh_h",deepcopy(mesh.points),mesh.lines,mesh.triangles,mesh.tetrahedra,mesh.domains,"constructed from mesh", mesh.tri2tet,[])
#     mesh_h.points[crdnt,pnt_idx]+=h
#     G=discretize(mesh_h, dscrp, C)
#     sol, n, flag = householder(G,ω0,maxiter=4, output = false, n_eig_val=3)
#     sens_FD[pnt_idx]=(sol.params[:ω]-ω0)/(h)
#
#     #update progressMeter
#     ProgressMeter.next!(p)
# end

## output
mode="[$(string(round((ω0/2/pi),digits=2)))]Hz"

data=Dict()
data["abs_p"*mode]=abs.(sol.v0)./maximum(abs.(sol.v0))
data["phase_p"*mode]=angle.(sol.v0)
data["abs_p_adj"*mode]=abs.(sol.v0Adj)./maximum(abs.(sol.v0Adj))
data["phase_p_adj"*mode]=angle.(sol.v0Adj)
data["discrete adjoint x"*mode]=real.(sens1)
data["discrete adjoint y"*mode]=real.(sens2)
data["discrete adjoint z"*mode]=real.(sens3)
#data["finite difference"*mode]=real.(sens_FD)
#data["central finite difference"]=real.(sens_CFD)
#data["diff DA -- FD"*mode]=abs.((sens.-sens_FD))
vtk_write("shape_sense_v1_fine", mesh, data)

## some tests
h=-1E-8
mesh2=deepcopy(mesh)
mesh2.points[3,surface_points[360]]+=h
L2=discretize(mesh2, dscrp, C)
##
sol2, n, flag = householder(L2,ω0,maxiter=4, output = false, n_eig_val=2)

##
sens1=-v0Adj'*D(sol.params[:ω])*v0

sens2=(sol2.params[:ω]-sol.params[:ω])/h

##    # #central finite difference failed because more elements need to be evaluated
    #sens_CFD=zeros(ComplexF64,size(mesh.points,2))
    # mesh_h=Mesh("mesh_h",deepcopy(mesh.points),mesh.lines,mesh.triangles,mesh.tetrahedra,domains,"constructed from mesh", mesh.tri2tet,[])
    # D_center=discretize(mesh_h, dscrp, C, mass_weighting=false)
    # D=deepcopy(L)
    # for (idx,term) in enumerate(D_center.terms)
    #     I,J,V=SparseArrays.findnz(term.coeff)
    #     for (ii,jj,vv) in zip(I,J,V)
    #         D.terms[idx].coeff[ii,jj]-=vv
    #     end
    # end
    # for (idx,term) in enumerate(D_right.terms)
    #     I,J,V=SparseArrays.findnz(term.coeff)
    #     for (ii,jj,vv) in zip(I,J,V)
    #         D.terms[idx].coeff[ii,jj]+=vv
    #     end
    # end
    # sol, n, flag = householder(D,ω0,maxiter=4, output = false, n_eig_val=3)
    # ω_right=sol.params[:ω]
    # D=deepcopy(L)
    # for (idx,term) in enumerate(D_center.terms)
    #     I,J,V=SparseArrays.findnz(term.coeff)
    #     for (ii,jj,vv) in zip(I,J,V)
    #         D.terms[idx].coeff[ii,jj]-=vv
    #     end
    # end
    # for (idx,term) in enumerate(D_left.terms)
    #     I,J,V=SparseArrays.findnz(term.coeff)
    #     for (ii,jj,vv) in zip(I,J,V)
    #         D.terms[idx].coeff[ii,jj]+=vv
    #     end
    # end
    # sol, n, flag = householder(D,ω0,maxiter=4, output = false, n_eig_val=3)
    # ω_left=sol.params[:ω]
    # sens_CFD[pnt_idx]=(ω_right-ω_left)/(2*h)

## post production
findmax(abs.(sens-sens_FD))
#(0.018001506831495817, CartesianIndex(2, 772))
# 0.00844 (1,8) #central difference h=default
#(0.4084884704070646, CartesianIndex(3, 553)) h=1E-6
# (6.740726362039766, CartesianIndex(2, 198)) h=1E-12
# (0.08539897971817823, CartesianIndex(1, 549)) h=1E-10
