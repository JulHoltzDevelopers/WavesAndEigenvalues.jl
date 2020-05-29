#module APE
using WavesAndEigenvalues.Meshutils, WavesAndEigenvalues.NLEVP
import SparseArrays, LinearAlgebra, ProgressMeter
include("./src/FEM.jl")

#
function potflow(mesh::Mesh,dscrp)
    N_points=size(mesh.points)[2]
    dim=deepcopy(N_points)
    MM=Float64[]
    II=Int64[]
    JJ=Int64[]

    for smplx in mesh.tetrahedra
        J=CooTrafo(mesh.points[:,smplx])
        ii, jj =create_indices(smplx)
        mm=s43nv1nu1(J)
        append!(MM,mm[:])
        append!(II,ii[:])
        append!(JJ,jj[:])
    end
    M=SparseArrays.sparse(II,JJ,MM,dim,dim)

    v=zeros(dim)
    for (dom,val) in dscrp
        compute_size!(mesh,dom)
        a=val/mesh.domains[dom]["size"]
        simplices=mesh.domains[dom]["simplices"]
        MM=Float64[]
        II=Int64[]
        for smplx in mesh.triangles[simplices]
            J=CooTrafo(mesh.points[:,smplx])
            mm=-s33v1(J)*a
            v[smplx[:]]+=mm[:]
        end
    end
    return M,v
end


##
function discretize(mesh::Mesh)
    L=LinearOperatorFamily(["ω","λ"],complex([0.,Inf]))
    N_points=size(mesh.points)[2]
    dim=deepcopy(N_points)
    P=101325
    ρ=1.225
    γ=1.4
    Y=1E15

    #Boundary matrix
    L.params[:Y]=Y
    MM=ComplexF64[]
    II=Int64[]
    JJ=Int64[]
    smplx_list=mesh.domains["Outlet"]["simplices"]
    for tri in mesh.triangles[smplx_list]
        J=CooTrafo(mesh.points[:,tri])
        ii, jj =create_indices(tri)
        mm=sqrt(γ*P/ρ)*s33v1u1(J)
        append!(MM,mm[:])
        append!(II,ii[:].+3*dim)
        append!(JJ,jj[:].+3*dim)
    end
        M=SparseArrays.sparse(II,JJ,MM,4*dim,4*dim)
    push!(L,Term(M,(pow1,),((:Y,),),"Y","B"))




    #Term II and IV
    MM=ComplexF64[]
    II=Int64[]
    JJ=Int64[]
    for tet in mesh.tetrahedra
        J=CooTrafo(mesh.points[:,tet])
        ii, jj =create_indices(tet)


        #TERM II
        #u†p
        mm=s43v1du1(J,1)
        append!(MM,mm[:])
        append!(II,ii[:].+0*dim)
        append!(JJ,jj[:].+3*dim)

        #v†p
        mm=s43v1du1(J,2)
        append!(MM,mm[:])
        append!(II,ii[:].+1*dim)
        append!(JJ,jj[:].+3*dim)

        #w†p
        mm=s43v1du1(J,3)
        append!(MM,mm[:])
        append!(II,ii[:].+2*dim)
        append!(JJ,jj[:].+3*dim)


        #termIV
        mm=-γ*P*s43dv1u1(J,1) #TODO: space-dependent P
        append!(MM,mm[:])
        append!(II,ii[:].+3*dim)
        append!(JJ,jj[:].+0*dim)
        mm=-γ*P*s43dv1u1(J,2) #TODO: space-dependent P
        append!(MM,mm[:])
        append!(II,ii[:].+3*dim)
        append!(JJ,jj[:].+1*dim)
        mm=-γ*P*s43dv1u1(J,3) #TODO: space-dependent P
        append!(MM,mm[:])
        append!(II,ii[:].+3*dim)
        append!(JJ,jj[:].+2*dim)
    end
    M=SparseArrays.sparse(II,JJ,MM,4*dim,4*dim)
    push!(L,Term(M,(),(),"","K"))

    #interior (term I + III)
    MM=ComplexF64[]
    II=Int64[]
    JJ=Int64[]
    #for tet in mesh.tetrahedra
    for tet in mesh.tetrahedra
        J=CooTrafo(mesh.points[:,tet])
        #velocity components term I
        mm=ρ*s43v1u1(J) #TODO: space-dependent ρ
        ii, jj =create_indices(tet)

        #u
        append!(MM,mm[:])
        append!(II,ii[:].+0*dim)
        append!(JJ,jj[:].+0*dim)
        #v
        append!(MM,mm[:])
        append!(II,ii[:].+1*dim)
        append!(JJ,jj[:].+1*dim)
        #w
        append!(MM,mm[:])
        append!(II,ii[:].+2*dim)
        append!(JJ,jj[:].+2*dim)

        #pressure component Term III
        mm=s43v1u1(J)
        append!(MM,mm[:])
        append!(II,ii[:].+3*dim)
        append!(JJ,jj[:].+3*dim)
    end

    M=SparseArrays.sparse(II,JJ,MM,4*dim,4*dim)
    push!(L,Term(1.0im*M,(pow1,),((:ω,),),"iω","M"))




    #grid mass matrix
    for tet in mesh.tetrahedra
        J=CooTrafo(mesh.points[:,tet])
        mm=s43v1u1(J)
        ii, jj =create_indices(tet)

        #u
        append!(MM,mm[:])
        append!(II,ii[:].+0*dim)
        append!(JJ,jj[:].+0*dim)
        #v
        append!(MM,mm[:])
        append!(II,ii[:].+1*dim)
        append!(JJ,jj[:].+1*dim)
        #w
        append!(MM,mm[:])
        append!(II,ii[:].+2*dim)
        append!(JJ,jj[:].+2*dim)

        #pressure component Term III
        append!(MM,mm[:])
        append!(II,ii[:].+3*dim)
        append!(JJ,jj[:].+3*dim)
    end
        M=SparseArrays.sparse(II,JJ,MM,4*dim,4*dim)
    push!(L,Term(-M,(pow1,),((:λ,),),"-λ","__aux__"))

    return L#II,JJ,MM,M


end

#end #module
##
mesh=Mesh("./examples/tutorials/Rijke_mm.msh",scale=0.001)
##
dscrp=Dict()
dscrp["Outlet"]=1.0
dscrp["Inlet"]=-1.0
L,v=potflow(mesh,dscrp)
phi=L\v
for smplx in mesh.tetrahedra
    J=CooTrafo(mesh.points[:,smplx])
    sum(J.inv\phi[smplx] #ZODO: Hier weitermachen
end
##
L=discretize(mesh)
##
sol,nn,flag=householder(L,500*2*pi,output=true,maxiter=30,tol=1E-10);
##
tri=mesh.triangles[1]
X=mesh.points[:,tri]
J=zeros(3,3)
J[:,1:2]=X[:,1:end-1].-X[:,end]
n=LinearAlgebra.cross(J[:,1],J[:,2])
n./=LinearAlgebra.norm(n)
J[:,3]=n
