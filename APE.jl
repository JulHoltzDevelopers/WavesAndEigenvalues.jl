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
        dim=4*N_points+length(mesh.triangles)+length(inner_triangles)
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
function discretize(mesh::Mesh,U)
    L=LinearOperatorFamily(["s","λ"],complex([0.,Inf]))
    N_points=size(mesh.points)[2]
    N_lines=length(mesh.lines)
    dim=N_points+3*(N_points+N_lines)
    P=101325
    ρ=1.225
    γ=1.4
    Y=0*1E15
    triangles,tetrahedra=aggregate_elements(mesh,:2)

    #interior (term I + III)
    MM=ComplexF64[]
    II=Int64[]
    JJ=Int64[]
    #for tet in mesh.tetrahedra
    for tet in tetrahedra
        J=CooTrafo(mesh.points[:,tet[1:4]])
        #velocity components term I
        mm=ρ*s43v2u2(J) #TODO: space-dependent ρ
        ii, jj =create_indices(tet)
        iip, jjp =create_indices(tet[1:4],tet[1:4])

        #u
        for dd=1:3
            append!(MM,mm[:])
            append!(II,ii[:].+N_points.+(dd-1)*(N_points+N_lines))
            append!(JJ,jj[:].+N_points.+(dd-1)*(N_points+N_lines))
        end

        #pressure component Term III
        mm=s43v1u1(J)
        append!(MM,mm[:])
        append!(II,iip[:])
        append!(JJ,jjp[:])
    end

    M=SparseArrays.sparse(II,JJ,MM,dim,dim)
    push!(L,Term(M,(pow1,),((:s,),),"s","M"))



    #Boundary matrix
    for dom in ("Inlet","Outlet")
        if dom=="Inlet"
            Ysym=:Y_in
        elseif dom =="Outlet"
            Ysym=:Y_out
        end
        L.params[Ysym]=Y

        MM=ComplexF64[]
        II=Int64[]
        JJ=Int64[]
        smplx_list=mesh.domains[dom]["simplices"]
        for tri in mesh.triangles[smplx_list]
            J=CooTrafo(mesh.points[:,tri])
            ii, jj =create_indices(tri)
            mm=sqrt(γ*P/ρ)*s33v1u1(J)
            append!(MM,mm[:])
            append!(II,ii[:])
            append!(JJ,jj[:])
        end
        M=SparseArrays.sparse(II,JJ,MM,dim,dim)
        push!(L,Term(M,(pow1,),((Ysym,),),string(Ysym),"B"))
    end



    #Term II and IV
    MM=ComplexF64[]
    II=Int64[]
    JJ=Int64[]
    for tet in tetrahedra
        J=CooTrafo(mesh.points[:,tet[1:4]])
        ii, jj = create_indices(tet[1:end],tet[1:4])
        iip,jjp= create_indices(tet[1:4],tet[1:end])
        for dd=1:3
            #TERM II
            #u†p
            mm=s43v2du1(J,dd)
            append!(MM,mm[:])
            append!(II,ii[:].+N_points.+(dd-1)*(N_points+N_lines))
            append!(JJ,jj[:])

            #termIV
            mm=-γ*P*s43dv1u2(J,dd) #TODO: space-dependent P
            append!(MM,mm[:])
            append!(II,iip[:])
            append!(JJ,jjp[:].+N_points.+(dd-1)*(N_points+N_lines))
        end
    end
    M=SparseArrays.sparse(II,JJ,MM,dim,dim)
    push!(L,Term(M,(),(),"","K"))




    #term V and VI (mean flow)
    MM=ComplexF64[]
    II=Int64[]
    JJ=Int64[]
    for tet in tetrahedra
        J=CooTrafo(mesh.points[:,tet[1:4]])
        ii, jj =create_indices(tet)
        iip, jjp =create_indices(tet[1:4])
        #term V
        for dd=1:3
            for ee=1:3
                u=U[ee,tet[1:4]]
                #mtx= #TODO: speed up by correct ordering
                mm=ρ*(s43diffc1(J,u,dd)*s43v2u2(J)+s43v2du2c1(J,u,dd))
                append!(MM,mm[:])
                append!(II,ii[:].+N_points.+(dd-1)*(N_points+N_lines))
                append!(JJ,jj[:].+N_points.+(ee-1)*(N_points+N_lines))
            end
            #term VI
            u=U[dd,tet[1:4]]
            mm=s43v1du1c1(J,u,dd)
            append!(MM,mm[:])
            append!(II,iip[:])
            append!(JJ,jjp[:])
        end
    end
    L.params[:v]=0.0
    M=SparseArrays.sparse(II,JJ,MM,dim,dim)
    push!(L,Term(M,(pow1,),((:v,),),"v","U"))



    #grid mass matrix
    MM=ComplexF64[]
    II=Int64[]
    JJ=Int64[]
    for tet in tetrahedra
        J=CooTrafo(mesh.points[:,tet[1:4]])
        mm=ρ*s43v2u2(J)
        ii, jj =create_indices(tet)
        iip, jjp =create_indices(tet[1:4])

        for dd=1:3
            #u
            append!(MM,mm[:])
            append!(II,ii[:].+N_points.+(dd-1)*(N_points+N_lines))
            append!(JJ,jj[:].+N_points.+(dd-1)*(N_points+N_lines))
        end
        #pressure component Term III
        mm=s43v1u1(J)
        append!(MM,mm[:])
        append!(II,iip[:])
        append!(JJ,jjp[:])
    end
    M=SparseArrays.sparse(II,JJ,MM,dim,dim)
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
L,v=potflow(mesh,dscrp,:h)
phi=L\v
##
N_points=size(mesh.points,2)
u=Array{Float64}(undef,3,N_points)
# u[1,:]=zeros(N_points)#phi[0*N_points+1:1*N_points]
# u[2,:]=zeros(N_points)#phi[1*N_points+1:2*N_points]
# u[3,:]=ones(N_points)#phi[2*N_points+1:3*N_points]
u[1,:]=phi[1*N_points+1:2*N_points]
u[2,:]=phi[2*N_points+1:3*N_points]
u[3,:]=phi[3*N_points+1:4*N_points]
A=mesh.domains["Inlet"]["size"]

##
L=discretize(mesh,u)
##
L.params[:Y_in]=0*1E15
L.params[:Y_out]=0*1E15
L.params[:v]=10
##
import Arpack
##
Y_in=0
Y_out=0
V=-.01
M=L.terms[1].coeff
A=Y_in*L.terms[2].coeff+Y_out*L.terms[3].coeff+L.terms[4].coeff+V*L.terms[5].coeff


SOL=Arpack.eigs(A,M,nev=5, sigma = 340im*2*pi,)

##
Γ=[335.0+5.0im, 335.0-5.0im, 345.0-5.0im, 345.0+5.0im].*2*pi*1im #corner points for the contour (in this case a rectangle)
Ω, P = beyn(L,Γ,l=5,N=16, output=true)
##
L.params[:v]=1
sol,nn,flag=householder(L,( 342.4168944643069im - 0.23084620589322044)*2*pi,output=true,maxiter=10,tol=1E-10);
##
N_points=size(mesh.points,2)
N_lines=length(mesh.lines)
data=Dict()
dim=size(mesh.points,2)
max_uvw=maximum(abs.(sol.v[N_points+1:end]))
p=sol.v[N_points.+(0)*(N_points+N_lines)+1:N_points.+(1)*(N_points+N_lines)]
data["abs u"]=abs.(p)/max_uvw
data["phase u"]=angle.(p)
data["real u"]=real.(p)
data["imag u"]=imag.(p)
p=sol.v[N_points.+(1)*(N_points+N_lines)+1:N_points.+(2)*(N_points+N_lines)]
data["abs v"]=abs.(p)/max_uvw
data["phase v"]=angle.(p)
data["real v"]=real.(p)
data["imag v"]=imag.(p)
p=sol.v[N_points.+(2)*(N_points+N_lines)+1:N_points.+(3)*(N_points+N_lines)]
data["abs w"]=abs.(p)/max_uvw
data["phase w"]=angle.(p)
data["real w"]=real.(p)
data["imag w"]=imag.(p)
p=sol.v[1:N_points]
data["abs p"]=abs.(p)/maximum(abs.(p))
data["phase p"]=angle.(p)
data["real p"]=real.(p)
data["imag p"]=imag.(p)
data["U"]=u[1,:]
data["V"]=u[2,:]
data["W"]=u[3,:]
vtk_write("flow",mesh,data)
##
perturb_fast!(sol,L,:v,30)
f(v)=sol(:v,v,15,15)

###
L.params[:v]=A*3
sol1,nn,flag=householder(L,370*2*pi,output=true,maxiter=10,tol=1E-10);
