module APE
#header
import SparseArrays, LinearAlgebra, ProgressMeter
using ..Meshutils, ..NLEVP
include("Meshutils_exports.jl")
include("NLEVP_exports.jl")
include("./FEM/FEM.jl")
export compute_potflow_field, discretize

function discretize(mesh::Mesh,dscrp,U;output=true)
    L=LinearOperatorFamily(["s","λ"],complex([0.,Inf]))
    N_points=size(mesh.points)[2]
    N_lines=length(mesh.lines)
    dim=N_points+3*(N_points+N_lines)
    #gas properties (air)
    P=101325 #ambient pressure (one atmosphere)
    ρ=1.225 #density
    γ=1.4 #ratio of specific heats

    triangles,tetrahedra=aggregate_elements(mesh,:quad)


    #initialize progress bar
    if output
        n_task=4*length(mesh.tetrahedra)
        for dom in keys(dscrp)
            n_task+=length(mesh.domains[dom]["simplices"])
        end

        prog = ProgressMeter.Progress(n_task,desc="Discretize... ", dt=.5,
             barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=20)

    end

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
        if output
            ProgressMeter.next!(prog)
        end
    end

    M=SparseArrays.sparse(II,JJ,MM,dim,dim)
    push!(L,Term(M,(pow1,),((:s,),),"s","M"))



    #Boundary matrix
    for (dom,val) in dscrp
        if dom=="Inlet"
            Ysym=:Y_in
        elseif dom =="Outlet"
            Ysym=:Y_out
        end
        L.params[Ysym]=-sqrt(γ*P/ρ)/(val/mesh.domains[dom]["size"])

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
            if output
                ProgressMeter.next!(prog)
            end
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
        if output
            ProgressMeter.next!(prog)
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
        if output
            ProgressMeter.next!(prog)
        end
    end
    L.params[:v]=1.0
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
        if output
            ProgressMeter.next!(prog)
        end
    end
    M=SparseArrays.sparse(II,JJ,MM,dim,dim)
    push!(L,Term(-M,(pow1,),((:λ,),),"-λ","__aux__"))
    return L#II,JJ,MM,M
end

"""
    U=compute_potflow_field(mesh::Mesh,dscrp;order=:lin,output=true)

Compute potential flow field `U` on mesh `mesh` using the boundary conditions
specified in `dscrp`.

# Arguments
- `mesh::Mesh`: discretization of the computational domain
- `dscrp::Dict`: Dictionairy mapping domain names to volume flow values. Positive and negative volume flow correspond to inflow and outflow, respectively.
- `order::Symbol=:lin`: (optional) toggle whether computed velocity field is constant on a tetrahedron (`:const`) or linear interpolated from the vertices (`:lin`).
- `output::Bool=false`: (optional) toggle showing progress bar.

# Returns
- `U::Array`: 3×`N`-Array. Each column corresponds to a velocity vector. Depending on whether "order==:const" or "order==:lin" the number of columns will be `N==size(mesh.points,2)` or `N==length(mesh.tetrahedra)`, respectively.

# Notes
The flow is computed by solving the Poisson equation. The specified volume flows
must sum up to 0, otherwise the continuum equation is violated.
"""
function compute_potflow_field(mesh::Mesh,dscrp;order=:lin,output=false)
    #the potential is a stem function of the velocity
    #the element order for the potential must therefore be one order higher
    if order==:const
        triangles, tetrahedra, dim = aggregate_elements(mesh,:lin)
    elseif order==:lin
        triangles, tetrahedra, dim = aggregate_elements(mesh,:herm)
    else
        println("Error: order $order not supported for potential flow!")
        return nothing #force crash
    end
    function s43nvnu(J)
        if order==:const
            return s43nv1nu1(J)
        #elseif order==:quad
        #    return s43nv2nu2(J)
        elseif order==:lin
            return s43nvhnuh(J)
        else
            return nothing #force crash
        end
    end
    function s33v(J)
        if order==:const
            return s33v1(J)
        #elseif order==:quad
        #    return s33v2(J)
    elseif order==:lin
            return s33vh(J)
        else
            return nothing #force crash
        end
    end


    MM=Float64[]
    II=Int64[]
    JJ=Int64[]

    if output
        n_task=length(mesh.tetrahedra)
        for dom in keys(dscrp)
            n_task+=length(mesh.domains[dom]["simplices"])
        end

        p = ProgressMeter.Progress(n_task,desc="Discretize... ", dt=.5,
             barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=20)

    end

    #discretize stiffness matrix (Laplacian)
    for smplx in tetrahedra
        J=CooTrafo(mesh.points[:,smplx[1:4]])
        ii, jj =create_indices(smplx)
        mm=s43nvnu(J)
        append!(MM,mm[:])
        append!(II,ii[:])
        append!(JJ,jj[:])
        if output
            ProgressMeter.next!(p)
        end
    end
    M=SparseArrays.sparse(II,JJ,MM,dim,dim)

    #discretize source terms from B.C.
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
            if output
                ProgressMeter.next!(p)
            end
        end
    end

    #solve for potential
    phi=M\v
    #compute velocity from grad phi
    if order==:const
        U=Array{Float64}(undef,3,length(mesh.tetrahedra))
        D=[1 0 0;
            0 1 0;
            0 0 1;
            -1 -1 -1]
        for (idx,tet) in enumerate(mesh.tetrahedra)
            J=CooTrafo(mesh.points[:,tet])
            U[:,idx]=transpose(phi[tet])*D*J.inv
        end
    elseif order==:lin
        #potential field was calculated with hermitian elements
        #both the potential and its derivative are part of the solution vector
        N_points=size(mesh.points,2)
        U=Array{Float64}(undef,3,N_points)
        U[1,:]=phi[1*N_points+1:2*N_points]
        U[2,:]=phi[2*N_points+1:3*N_points]
        U[3,:]=phi[3*N_points+1:4*N_points]
    end
    return  U
end

end
