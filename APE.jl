#module APE
using WavesAndEigenvalues.Meshutils, WavesAndEigenvalues.NLEVP
import WavesAndEigenvalues.Meshutils: find_smplx, insert_smplx!
import SparseArrays, LinearAlgebra, ProgressMeter
import WavesAndEigenvalues
include("./src/FEM.jl")


function collect_triangles(mesh::Mesh)
    inner_triangles=[]
    for tet in mesh.tetrahedra
        for tri in [tet[[1,2,3]],tet[[1,2,4]],tet[[1,3,4]],tet[[2,3,4]] ]
            idx=find_smplx(mesh.triangles,tri)
            if idx==0
                insert_smplx!(inner_triangles,tri)
            end
        end
    end
    return inner_triangles
end

##

##


import ..Meshutils: get_line_idx
function aggregate_elements(mesh::Mesh, el_type=:1)
    N_points=size(mesh.points)[2]
    if (el_type in (:2,:h) ) &&  length(mesh.lines)==0
        collect_lines!(mesh)
    end

    if el_type==:1
        tetrahedra=mesh.tetrahedra
        triangles=mesh.triangles
        dim=N_points
    elseif el_type==:2
        triangles=Array{Array{UInt32,1}}(undef,length(mesh.triangles))
        tetrahedra=Array{Array{UInt32,1}}(undef,length(mesh.tetrahedra))
        tet=Array{UInt32}(undef,10)
        tri=Array{UInt32}(undef,6)
        for (idx,smplx) in enumerate(mesh.tetrahedra)
            tet[1:4]=smplx[:]
            tet[5]=get_line_idx(mesh,smplx[[1,2]])+N_points#find_smplx(mesh.lines,smplx[[1,2]])+N_points #TODO: type stability
            tet[6]=get_line_idx(mesh,smplx[[1,3]])+N_points
            tet[7]=get_line_idx(mesh,smplx[[1,4]])+N_points
            tet[8]=get_line_idx(mesh,smplx[[2,3]])+N_points
            tet[9]=get_line_idx(mesh,smplx[[2,4]])+N_points
            tet[10]=get_line_idx(mesh,smplx[[3,4]])+N_points
            tetrahedra[idx]=copy(tet)
        end
        for (idx,smplx) in enumerate(mesh.triangles)
            tri[1:3]=smplx[:]
            tri[4]=get_line_idx(mesh,smplx[[1,2]])+N_points
            tri[5]=get_line_idx(mesh,smplx[[1,3]])+N_points
            tri[6]=get_line_idx(mesh,smplx[[2,3]])+N_points
            triangles[idx]=copy(tri)
        end
        dim=N_points+length(mesh.lines)
    elseif el_type==:h
        #if mesh.tri2tet[1]==0xffffffff
        #    link_triangles_to_tetrahedra!(mesh)
        #end
        inner_triangles=collect_triangles(mesh) #TODO integrate into mesh structure
        triangles=Array{Array{UInt32,1}}(undef,length(mesh.triangles))
        tetrahedra=Array{Array{UInt32,1}}(undef,length(mesh.tetrahedra))
        tet=Array{UInt32}(undef,20)
        tri=Array{UInt32}(undef,13)
        for (idx,smplx) in enumerate(mesh.triangles)
            tri[1:3]  =  smplx[:]
            tri[4:6]  =  smplx[:].+N_points
            tri[7:9]  =  smplx[:].+2*N_points
            tri[10:12]  =  smplx[:].+3*N_points
            fcidx     =  find_smplx(mesh.triangles,smplx)
            if fcidx !=0
                tri[13]   =  fcidx+4*N_points
            else
                tri[13]   =  find_smplx(inner_triangles,smplx)+4*N_points+length(mesh.triangles)
            end
            triangles[idx]=copy(tri)
        end
        for (idx, smplx) in enumerate(mesh.tetrahedra)
            tet[1:4] = smplx[:]
            tet[5:8] = smplx[:].+N_points
            tet[9:12] = smplx[:].+2*N_points
            tet[13:16]= smplx[:].+3*N_points

            for (jdx,tria) in enumerate([smplx[[2,3,4]],smplx[[1,3,4]],smplx[[1,2,4]],smplx[[1,2,3]]])
                fcidx     =  find_smplx(mesh.triangles,tria)
                if fcidx !=0
                    tet[16+jdx]   =  fcidx+4*N_points
                else
                    fcidx=find_smplx(inner_triangles,tria)
                    if fcidx==0
                        println("Error, face not found!!!")
                        return nothing
                    end
                    tet[16+jdx]   =  fcidx+4*N_points+length(mesh.triangles)
                end
            end
            tetrahedra[idx]=copy(tet)
        end
        dim=N_points*4+length(mesh.triangles)+length(inner_triangles)
    end
    return triangles, tetrahedra, dim
end


function potflow(mesh::Mesh,dscrp,order)
    triangles, tetrahedra, dim = aggregate_elements(mesh,order)
    function s43nvnu(J)
        if order==:1
            return s43nv1nu1(J)
        elseif order==:2
            return s43nv2nu2(J)
        elseif order==:h
            return s43nvhnuh(J)
        else
            return nothing #force crash
        end
    end
    function s33v(J)
        if order==:1
            return s33v1(J)
        elseif order==:2
            return s33v2(J)
        elseif order==:h
            return s33vh(J)
        else
            return nothing #force crash
        end
    end


    MM=Float64[]
    II=Int64[]
    JJ=Int64[]

    for smplx in tetrahedra
        J=CooTrafo(mesh.points[:,smplx[1:4]])
        ii, jj =create_indices(smplx)
        mm=s43nvnu(J)
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
        for smplx in triangles[simplices]
            J=CooTrafo(mesh.points[:,smplx[1:3]])
            mm=-s33v(J)*a
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
    Y=0*1E15

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
L,v=potflow(mesh,dscrp,:1)
phi=L\v
##
data=Dict()
data["p1"]=phi1[1:size(mesh.points,2)]
data["p2"]=phi2[1:size(mesh.points,2)]
data["ph"]=phih[1:size(mesh.points,2)]
vtk_write("potential",mesh,data)

##
for smplx in mesh.tetrahedra
    J=CooTrafo(mesh.points[:,smplx])
    sum(J.inv\phi[smplx] #ZODO: Hier weitermachen
end
##
L=discretize(mesh)
##
sol,nn,flag=householder(L,670*2*pi,output=true,maxiter=30,tol=1E-10);
##
tri=mesh.triangles[1]
X=mesh.points[:,tri]
J=zeros(3,3)
J[:,1:2]=X[:,1:end-1].-X[:,end]
n=LinearAlgebra.cross(J[:,1],J[:,2])
n./=LinearAlgebra.norm(n)
J[:,3]=n
