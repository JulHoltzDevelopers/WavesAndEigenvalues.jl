"""
Module providing functionality to numerically discretize the (thermoacoustic) Helmholtz equation by first and second-order finite elements.
"""
module Helmholtz
import SparseArrays
using ..Meshutils, ..NLEVP
include("Meshutils_exports.jl")
include("NLEVP_exports.jl")
export discretize
##
function outer(ii,aa,jj,bb)
    #TODO: preallocation
    mm=ComplexF64[]
    II=UInt32[]
    JJ=UInt32[]
    for (a,i) in zip(aa,ii)
        for (b,j) in zip(bb,jj)
            push!(mm,a*b)
            push!(II,i)
            push!(JJ,j)
        end
    end
    #TODO: conversion to array
    return II,JJ,mm
end



function blochify(ii,jj,mm, naxis,nsector; axis = true)
    #nsector=naxis+nxbloch+nbody+nxsymmetry+nbody
    blochshift=nsector-naxis
    MM=ComplexF64[]
    MM_plus=ComplexF64[]
    MM_minus=ComplexF64[]
    II=UInt32[]
    JJ=UInt32[]
    II_plus=UInt32[]#TODO: consider preallocation
    JJ_plus=UInt32[]
    II_minus=UInt32[]
    JJ_minus=UInt32[]
    #println("$(length(ii)) $(length(jj)) $(length(mm))")
    for (i,j,m) in zip(ii,jj,mm)
        #TODO: Hier eine if abfrage ob dof <= npoints dann kann man die lines gesondert behandeln.
        i_check = i>nsector
        j_check = j>nsector

        # map Bloch_ref to Bloch_img
        if i_check
            i-=blochshift
        end
        if j_check
            j-=blochshift
        end

        if axis &&( i<=naxis || j<=naxis) &&  !(i<=naxis && j<=naxis)
            m=0
        end

        if (!i_check && !j_check) || (i_check && j_check) #no manipulation of matrix entries
            append!(II,i)
            append!(JJ,j)
            append!(MM,m)

        elseif !i_check && j_check
            append!(II_plus,i)
            append!(JJ_plus,j)
            append!(MM_plus,m)

        elseif i_check && !j_check
            append!(II_minus,i)
            append!(JJ_minus,j)
            append!(MM_minus,m)
        else
            println("ERROR: in blochification")
            return
        end
    end

    #println("Here:::",length(II)," ", length(JJ)," ", length(MM))
    return (II,II_plus,II_minus), (JJ,JJ_plus,JJ_minus), (MM, MM_plus, MM_minus)
end


##
module lin
export assemble_volume_source, assemble_gradient_source, assemble_mass_operator, assemble_stiffness_operator, assemble_boundary_mass_operator
include("FEMlin.jl")
end
module quad20
export assemble_volume_source, assemble_gradient_source, assemble_mass_operator, assemble_stiffness_operator, assemble_boundary_mass_operator
include("FEMquadlin.jl")
end

module quad
export assemble_volume_source, assemble_gradient_source, assemble_mass_operator, assemble_stiffness_operator, assemble_boundary_mass_operator
include("FEMquad.jl")
end

import .lin
import .quad
import .quad20

"""
    L=discretize(mesh, dscrp, C; el_type=1, c_type=0, b=:__none__)

Discretize the Helmholtz equation using the mesh `mesh`.

# Arguments
- `mesh::Mesh`: tetrahedral mesh
- `dscrp::Dict `: dictionary containing information on where to apply physical constraints. Such as standard wave propagation, boundary conditions, flame responses, etc.
- `C:Array`: array defining the speed of sound. If "c_type==0" the speed of sound is constant along one tetrahedron and `length(C)==length(mesh.tetrahedra)`. If `c_type==1` the speed of sound is linearly interpolated between the vertices of the mesh and `length(C)==size(mesh.points,2)`.
- `eltype = 1`: optional paramater to select between first (`eltype==1` the default) and second order (`eltype==2`) finite elements.
- `c_type = 1`: optional parameter controlling whether speed of sound is constant on a tetrahedron or linearly interpolated between vertices.
- `b::Symbol=:__none__`: optional parameter defining the Bloch wave number. If `b=:__none__` (the default) noch Blochwave formalism is applied.

# Returns
- `L::LinearOperatorFamily`: parametereized discretization of the specified Helmholtz equation.
"""
function discretize(mesh::Mesh, dscrp, C; el_type=1, c_type=0, b=:__none__)
    if el_type==1 && c_type==0
        println("using linear finite elements and locally constant speed of sound")
        assemble_volume_source, assemble_gradient_source, assemble_mass_operator, assemble_stiffness_operator, assemble_boundary_mass_operator =lin.assemble_volume_source, lin.assemble_gradient_source, lin.assemble_mass_operator, lin.assemble_stiffness_operator, lin.assemble_boundary_mass_operator
    elseif el_type==2 && c_type==0
        println("using linear finite elements and locally constant speed of sound")
        assemble_volume_source, assemble_gradient_source, assemble_mass_operator, assemble_stiffness_operator, assemble_boundary_mass_operator =quad20.assemble_volume_source, quad20.assemble_gradient_source, quad20.assemble_mass_operator, quad20.assemble_stiffness_operator, quad20.assemble_boundary_mass_operator
    elseif el_type==2 && c_type==1
        println("using quadratic finite elements and locally linear speed of sound")
                assemble_volume_source, assemble_gradient_source, assemble_mass_operator, assemble_stiffness_operator, assemble_boundary_mass_operator =quad.assemble_volume_source, quad.assemble_gradient_source, quad.assemble_mass_operator, quad.assemble_stiffness_operator, quad.assemble_boundary_mass_operator
    else
        println("ERROR: Elements not supported!")
        return
    end


    L=LinearOperatorFamily(["ω","λ"],complex([0.,Inf]))
    dim=size(mesh.points)[2]
    if b!=:__none__
        bloch=true
        naxis=mesh.dos[2]
        nsector=naxis+mesh.dos[3]+2*mesh.dos[4]+mesh.dos[5]#naxis+nxbloch+nbody+nxsymmetry+nbody
        Δϕ=2*pi/mesh.dos[1]
        exp_plus(z,k)=exp_az(z,Δϕ*1.0im,k) #check whether this is a closure
        exp_minus(z,k)=exp_az(z,-Δϕ*1.0im,k)
        txt_plus="*exp(i$(b)2π/$(mesh.dos[1]))"
        txt_minus="*exp(-i$(b)2π/$(mesh.dos[1]))"
        L.params[b]=0.0+0.0im
        dim-=mesh.dos[3]
    else
        bloch=false
        naxis=nsector=0
    end


    #aggregator to create local elements ( tetrahedra, triangles)
    if el_type==2
        N_points=UInt32(dim)
        dim+=size(mesh.lines)[1] #TODO: check for bloch
        triangles=Array{UInt32,1}[] #TODO: preallocation?
        tetrahedra=Array{UInt32,1}[]
        tet=Array{UInt32}(undef,10)
        tri=Array{UInt32}(undef,6)
        for (idx,smplx) in enumerate(mesh.tetrahedra) #TODO: no enumeration
            tet[1:4]=smplx[:]
            tet[5]=find_smplx(mesh.lines,smplx[[1,2]])+N_points #TODO: type stability
            tet[6]=find_smplx(mesh.lines,smplx[[1,3]])+N_points
            tet[7]=find_smplx(mesh.lines,smplx[[1,4]])+N_points
            tet[8]=find_smplx(mesh.lines,smplx[[2,3]])+N_points
            tet[9]=find_smplx(mesh.lines,smplx[[2,4]])+N_points
            tet[10]=find_smplx(mesh.lines,smplx[[3,4]])+N_points
            push!(tetrahedra,copy(tet))
        end
        for (idx,smplx) in enumerate(mesh.triangles)
            tri[1:3]=smplx[:]
            tri[4]=find_smplx(mesh.lines,smplx[[1,2]])+N_points
            tri[5]=find_smplx(mesh.lines,smplx[[1,3]])+N_points
            tri[6]=find_smplx(mesh.lines,smplx[[2,3]])+N_points
            push!(triangles,copy(tri))
        end

    elseif el_type==1
            tetrahedra=mesh.tetrahedra
            triangles=mesh.triangles
    end

    if c_type==0
        C_tet=C
        link_triangles_to_tetrahedra!(mesh)
        C_tri=C[mesh.tri2tet]
    elseif c_type==1
        C_tet=Array{UInt32,1}[] #TODO: Preallocation
        C_tri=Array{UInt32,1}[]
        for tet in tetrahedra
            push!(C_tet,C[tet[1:4]])
        end
        for tri in triangles
            push!(C_tri,C[tri[1:3]])
        end
    end

    ## build matrices from domain definitions
    for (domain,(type,data)) in dscrp
        if type==:interior
            make=[:M,:K]

        elseif type==:admittance
            make=[:C]
            adm_sym,adm_val=data
            adm_txt="ω*"*string(adm_sym)
            if adm_sym ∉ keys(L.params)
                L.params[adm_sym]=adm_val
            end


        elseif type==:flame
            make=[:Q]
            gamma,rho,nglobal,x_ref,n_ref,n_sym,tau_sym,n_val,tau_val=data
            nlocal=(gamma-1)/rho*nglobal/compute_volume!(mesh,domain)
            if n_sym ∉ keys(L.params)
                L.params[n_sym]=n_val
            end
            if tau_sym ∉ keys(L.params)
                L.params[tau_sym]=tau_val
            end
            flame_func=(pow1,exp_delay,)
            flame_arg=((n_sym,),(:ω,tau_sym))
            flame_txt="$(string(n_sym))*exp(-iω$(string(tau_sym)))"

        elseif type==:flameresponse
            make=[:Q]
            gamma,rho,nglobal,x_ref,n_ref,eps_sym,eps_val=data
            nlocal=(gamma-1)/rho*nglobal/compute_volume!(mesh,domain)
            if eps_sym ∉ keys(L.params)
                L.params[eps_sym]=eps_val
            end
            flame_func=(pow1,)
            flame_arg=((eps_sym,),)
            flame_txt="$(string(eps_sym))"
        elseif type==:fancyflame
            make=[:Q]
            gamma,rho,nglobal,x_ref,n_ref,n_sym,tau_sym, a_sym, n_val, tau_val, a_val=data
            nlocal=(gamma-1)/rho*nglobal/compute_volume!(mesh,domain)
            if n_sym ∉ keys(L.params)
                L.params[n_sym]=n_val
            end
            if tau_sym ∉ keys(L.params)
                L.params[tau_sym]=tau_val
            end
            if a_sym ∉ keys(L.params)
                L.params[a_sym]=a_val
            end

            flame_func=(pow1,exp_ax2mxit,)
            flame_arg=((n_sym,),(:ω,tau_sym,a_sym,))
            flame_txt="$(string(n_sym))* exp($(string(a_sym))ω^2-iω$(string(tau_sym)))"
        else
            make=[]
        end

        for opr in make
            if opr==:M
                I,J,V=assemble_mass_operator(mesh.points,tetrahedra[mesh.domains[domain]["simplices"]])
                func=(pow2,)
                arg=((:ω,),)
                txt="ω^2"
                mat="M"
            elseif opr==:K
                I,J,V=assemble_stiffness_operator(mesh.points,tetrahedra[mesh.domains[domain]["simplices"]],C_tet[mesh.domains[domain]["simplices"]])
                func=()
                arg=()
                txt=""
                #V=-V
                mat="K"
            elseif opr==:C
                I,J,V=assemble_boundary_mass_operator(mesh.points,triangles[mesh.domains[domain]["simplices"]],C_tri[mesh.domains[domain]["simplices"]])
                V*=-1im
                func=(pow1,pow1,)
                arg=((:ω,),(adm_sym,),)
                txt=adm_txt
                mat="C"
            elseif opr==:Q
                I,S=assemble_volume_source(mesh.points,tetrahedra[mesh.domains[domain]["simplices"]])
                ref_idx=find_tetrahedron_containing_point(mesh,x_ref)
                J,G=assemble_gradient_source(mesh.points,tetrahedra[ref_idx],x_ref,n_ref)
                G=-nlocal*G
                #G=-G
                I,J,V=outer(I,S,J,G)
                func =flame_func
                arg=flame_arg
                txt=flame_txt
                mat="Q"
            end

            if bloch
                for (i,j,v,f,a,t) in zip(blochify(I,J,V,naxis,nsector)...,[(),(exp_plus,),(exp_minus,)],[(), ((b,),), ((b,),)], ["", txt_plus, txt_minus])
                    M =SparseArrays.sparse(i,j,v,dim,dim)
                    push!(L,Term(M,(func...,f...),(arg...,a...),txt*t,mat))
                end

            else
                M =SparseArrays.sparse(I,J,V,dim,dim)
                push!(L,Term(M,func,arg,txt,mat))
            end

        end
    end

    I,J,V=assemble_mass_operator(mesh.points,tetrahedra)
    if bloch
        IJV=blochify(I,J,V,naxis,nsector,axis=false)#TODO check whether this matrix needs more blochfication... especially because tetrahedra is not periodic
        I=[IJV[1][1]..., IJV[1][2]..., IJV[1][3]...]
        J=[IJV[2][1]..., IJV[2][2]..., IJV[2][3]...]
        V=[IJV[3][1]..., IJV[3][2]..., IJV[3][3]...]
    end
    M=SparseArrays.sparse(I,J,-V,dim,dim)
    push!(L,Term(M,(pow1,),((:λ,),),"-λ","__aux__"))
    return L
end
end#module
