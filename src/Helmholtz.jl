"""
Module providing functionality to numerically discretize the (thermoacoustic) Helmholtz equation by first, second, and hermitian-order finite elements.
"""
module Helmholtz
import SparseArrays, LinearAlgebra, ProgressMeter
import FFTW: fft
using ..Meshutils, ..NLEVP #TODO:check where find_smplx is introduced to scope
import ..Meshutils: get_line_idx
import ..NLEVP: generate_1_gz
include("Meshutils_exports.jl")
include("NLEVP_exports.jl")
include("./FEM/FEM.jl")
include("shape_sensitivity.jl")
include("Bloch.jl")
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

"""
    L=discretize(mesh, dscrp, C; order=:lin, b=:__none__, mass_weighting=true,source=false, output=true)

Discretize the Helmholtz equation using the mesh `mesh`.

# Arguments
- `mesh::Mesh`: tetrahedral mesh
- `dscrp::Dict `: dictionary containing information on where to apply physical constraints. Such as standard wave propagation, boundary conditions, flame responses, etc.
- `C:Array`: array defining the speed of sound. If `length(C)==length(mesh.tetrahedra)` the speed of sound is constant along one tetrahedron. If `length(C)==size(mesh.points,2)` the speed of sound is linearly interpolated between the vertices of the mesh.
- `order::Symbol = :lin`: optional paramater to select between first (`order==:lin` the default), second (`order==:quad`),or hermitian-order (`order==:herm`) finite elements.
- `b::Symbol=:__none__`: optional parameter defining the Bloch wave number. If `b=:__none__` (the default) no Blochwave formalism is applied.
- `mass_weighting=true`: optional parameter if true mass matrix is used as weighting matrix for householder, otherwise this matrix is not set.
- `source::Bool=false`: optional parameter to toggle the return of a source vector (experimental)
- `output::Bool=false': optional parameter to toggle live progress report.

# Returns
- `L::LinearOperatorFamily`: parametereized discretization of the specified Helmholtz equation.
- `rhs::LinearOperatorFamily`: parameterized discretization of the source vector. Only returned if `source==true`.  (experimental)
"""
function discretize(mesh::Mesh, dscrp, C; order=:lin, b=:__none__, mass_weighting=true, source=false, output=true)
    triangles,tetrahedra,dim=aggregate_elements(mesh,order)
    N_points=size(mesh.points,2)


    if length(C)==length(mesh.tetrahedra)
        C_tet=C
        if length(mesh.tri2tet)!=0 && mesh.tri2tet[1]==0xffffffff #TODO: Initialize as empty instead with sentinel value
            link_triangles_to_tetrahedra!(mesh)
        end
        C_tri=C[mesh.tri2tet]
    elseif length(C)==size(mesh.points,2)
        C_tet=Array{Float64,1}[] #TODO: Preallocation
        C_tri=Array{Float64,1}[]
        for tet in tetrahedra
            push!(C_tet,C[tet[1:4]])
        end
        for tri in triangles
            push!(C_tri,C[tri[1:3]])
        end
    end


    #initialize linear operator family...
    L=LinearOperatorFamily(["ω","λ"],complex([0.,Inf]))
    #...and source vector
    rhs=LinearOperatorFamily(["ω"],complex([0.,]))

    if b!=:__none__
        bloch=true
        naxis=mesh.dos.naxis
        nxbloch=mesh.dos.nxbloch
        nsector=naxis+mesh.dos.nxsector
        naxis_ln=mesh.dos.naxis_ln+N_points
        nsector_ln=mesh.dos.naxis_ln+mesh.dos.nxsector_ln+N_points
        Δϕ=2*pi/mesh.dos.DOS
        exp_plus(z,k)=exp_az(z,Δϕ*1.0im,k) #check whether this is a closure
        exp_minus(z,k)=exp_az(z,-Δϕ*1.0im,k)
        bloch_filt=zeros(ComplexF64,mesh.dos.DOS)
        bloch_filt[1]=(1.0+0.0im)/mesh.dos.DOS
        bloch_filt=fft(bloch_filt)
        bloch_filt=generate_Σy_exp_ikx(bloch_filt)
        anti_bloch_filt=generate_1_gz(bloch_filt)
        bloch_exp_plus=generate_gz_hz(bloch_filt,exp_plus)
        bloch_exp_minus=generate_gz_hz(bloch_filt,exp_minus)
        txt_plus="*exp(i$(b)2π/$(mesh.dos.DOS))"
        txt_minus="*exp(-i$(b)2π/$(mesh.dos.DOS))"
        txt_filt="*δ($b)"
        txt_filt_plus=txt_filt*txt_plus
        txt_filt_minus=txt_filt*txt_minus
        #txt_filt="*1($b)"
        #bloch_filt=pow0
        L.params[b]=0.0+0.0im
        if order==:lin
            dim-=mesh.dos.nxbloch
        elseif order==:quad
            dim-=mesh.dos.nxbloch+mesh.dos.nxbloch_ln
        elseif order==:herm
            dim-=4*mesh.dos.nxbloch #TODO: check blochify for hermitian elements
        end
        #N_points_dof=N_points-mesh.dos.nxbloch
    else
        bloch=false
        naxis=nsector=0
    end
    #wrapper for FEM constructors TODO: use multipledispatch and move to FEM.jl
    function stiff(J,c)
        if order==:lin
            if length(c)==1
                return -c^2*s43nv1nu1(J)
            elseif length(c)==4
                return -s43nv1nu1cc1(J,c)
            end
        elseif order==:quad
            if length(c)==1
                return -c^2*s43nv2nu2(J)
            elseif length(c)==4
                return -s43nv2nu2cc1(J,c)
            end
        elseif order==:herm
            if length(c)==1
                return -c^2*s43nvhnuh(J)
            elseif length(c)==4
                return -s43nvhnuhcc1(J,c)
            end
        end
    end
    function mass(J)
        if order==:lin
            return s43v1u1(J)
        elseif order==:quad
            return s43v2u2(J)
        elseif order==:herm
            return s43vhuh(J)#massh(J)
        end
    end

    function bound(J,c)
        if order==:lin
            if length(c)==1
                return c*s33v1u1(J)
            elseif length(c)==3
                return s33v1u1c1(J,c)
            end
        elseif order==:quad
            if length(c)==1
                return c*s33v2u2(J)
            elseif length(c)==3
                return s33v2u2c1(J,c)
            end
        elseif order==:herm
            if length(c)==1
                return c*s33vhuh(J)
            elseif length(c)==3
                return s33vhuhc1(J,c)
            end
        end
    end

    function volsrc(J)
        if order==:lin
            return s43v1(J)
        elseif order==:quad
            return s43v2(J)
        elseif order==:herm
            return s43vh(J)
        end
    end

    function gradsrc(J,n,x)
        if order==:lin
            return s43nv1rx(J,n,x)
        elseif order==:quad
            return s43nv2rx(J,n,x)
        elseif order==:herm
            return s43nvhrx(J,n,x)
        end
    end

    function wallsrc(J,c)
        if length(c)==1
            if order==:lin
                return c*s33v1(J)
            elseif order==:quad
                return c*s33v2(J)
            elseif order==:herm
                return c*s33vh(J)
            end
        else
            if order==:lin
                return s33v1c1(J,c)
            elseif order==:quad
                return s33v2c1(J,c)
            elseif order==:herm
                return s33vhc1(J,c)
            end
        end
    end

    #prepare progressbar
    if output
        n_task=0
        for (domain,(type,data)) in dscrp
            if type==:interior
                task_factor=2
            else
                task_factor=1
            end
                n_task+=task_factor*length(mesh.domains[domain]["simplices"])
        end

        p = ProgressMeter.Progress(n_task,desc="Discretize... ", dt=.5,
             barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=20)

    end

    ## build discretization matrices and vectors from domain definitions
    for (domain,(type,data)) in dscrp
        if type==:interior
            make=[:M,:K]
            stiff_func=()
            stiff_arg=()
            stiff_txt=""

        elseif type==:mass
            make=[:M]

        elseif type==:stiff
            make=[:K]
            stiff_func,stiff_arg,stiff_txt=data
            for args in stiff_arg
                for arg in args
                    L.params[arg]=0.0
                end
            end

        elseif type in (:admittance,:speaker)
            make=[]
            if type==:speaker
                append!(make,[:m])
                speak_sym,speak_val = data[1:2]
                rhs.params[speak_sym]=speak_val
                data = data[3:end]
            end


            if length(data)>0
                append!(make,[:C])
                if length(data)==2
                    adm_sym,adm_val=data

                    adm_txt="ω*"*string(adm_sym)
                    if adm_sym ∉ keys(L.params)
                        L.params[adm_sym]=adm_val
                        if type==:speaker
                            rhs.params[adm_sym]=adm_val
                        end
                    end
                    boundary_func=(pow1,pow1,)
                    boundary_arg=((:ω,),(adm_sym,),)
                    boundary_txt=adm_txt
                elseif length(data)==1
                    boundary_func=(generate_z_g_z(data[1]),)
                    boundary_arg=((:ω,),)
                    boundary_txt="ω*Y(ω)"
                elseif length(data)==4
                    Ass,Bss,Css,Dss = data
                    adm_txt="ω*C_s(iωI-A)^{-1}B"
                    func_z_stsp = generate_z_g_z(generate_stsp_z(Ass,Bss,Css,Dss))
                    boundary_func=(func_z_stsp,)
                    boundary_arg=((:ω,),)
                    boundary_txt=adm_txt
                end
            end



        elseif type==:flame #TODO: unified interface with flameresponse and normalize n_ref
            make=[:Q]
            isntau=false
            if length(data)==9 ## ref_idx is not specified by the user
                gamma,rho,nglobal,x_ref,n_ref,n_sym,tau_sym,n_val,tau_val=data
                ref_idx=0
                isntau=true
            elseif length(data)==10## ref_idx is specified by the user
                gamma,rho,nglobal,ref_idx,x_ref,n_ref,n_sym,tau_sym,n_val,tau_val=data
                isntau=true
            elseif length(data)==6 #custom FTF
                gamma,rho,nglobal,x_ref,n_ref,FTF=data
                ref_idx=0
                flame_func = (FTF,)
                flame_arg = ((:ω,),)
                if applicable(FTF,:ω)
                    flame_txt=FTF(:ω)
                else
                    flame_txt = "FTF(ω)"
                end
            elseif length(data)==5 #plain FTF
                gamma,rho,nglobal,x_ref,n_ref=data
                ref_idx=0
                L.params[:FTF]=0.0
                flame_func = (pow1,)
                flame_arg = ((:FTF,),)
                flame_txt = "FTF"
                gamma,rho,nglobal,x_ref,n_ref=data
            else
             println("Error: Data length does not match :flame option!")
            end


            nlocal=(gamma-1)/rho*nglobal/compute_size!(mesh,domain)
            if isntau
                if n_sym ∉ keys(L.params)
                    L.params[n_sym]=n_val
                end
                if tau_sym ∉ keys(L.params)
                    L.params[tau_sym]=tau_val
                end
                flame_func=(pow1,exp_delay,)
                flame_arg=((n_sym,),(:ω,tau_sym))
                flame_txt="$(string(n_sym))*exp(-iω$(string(tau_sym)))"
            end

            if ref_idx==0
                ref_idx=find_tetrahedron_containing_point(mesh,x_ref)
            end
            if ref_idx ∈ mesh.domains[domain]["simplices"]
                println("Warning: your reference point is inside the domain of heat release. (short-circuited FTF!)")
            end

        elseif type==:flameresponse
            make=[:Q]
            gamma,rho,nglobal,x_ref,n_ref,eps_sym,eps_val=data
            ref_idx=0
            nlocal=(gamma-1)/rho*nglobal/compute_size!(mesh,domain)
            if eps_sym ∉ keys(L.params)
                L.params[eps_sym]=eps_val
            end
            flame_func=(pow1,)
            flame_arg=((eps_sym,),)
            flame_txt="$(string(eps_sym))"
            if ref_idx==0
                ref_idx=find_tetrahedron_containing_point(mesh,x_ref)
            end



        elseif type==:fancyflame
            make=[:Q]
            gamma,rho,nglobal,x_ref,n_ref,n_sym,tau_sym, a_sym, n_val, tau_val, a_val=data
            ref_idx=0
            nlocal=(gamma-1)/rho*nglobal/compute_size!(mesh,domain)

            if typeof(n_val)<:Number #TODO: sanity check that other lists are same length
                if n_sym ∉ keys(L.params)
                    L.params[n_sym]=n_val
                end
                if tau_sym ∉ keys(L.params)
                    L.params[tau_sym]=tau_val
                end
                if a_sym ∉ keys(L.params)
                    L.params[a_sym]=a_val
                end

                flame_func=(pow1,exp_az2mzit,)
                flame_arg=((n_sym,),(:ω,tau_sym,a_sym,))
                flame_txt="$(string(n_sym))* exp($(string(a_sym))ω^2-iω$(string(tau_sym)))"
            else
                flame_arg=(:ω,)
                flame_txt=""
                for (ns, ts, as, nv, tv, av) in zip(n_sym, tau_sym, a_sym, n_val, tau_val,a_val)
                    L.params[ns]=nv
                    L.params[ts]=tv
                    L.params[as]=av
                    flame_arg=(flame_arg...,ns,ts,as,)
                    flame_txt*="[$(string(ns))* exp($(string(as))ω^2-iω$(string(ts)))+"
                end
                flame_txt=flame_txt[1:end-1]*"]"
                flame_arg=(flame_arg,)
                flame_func=(Σnexp_az2mzit,)
            end
            if ref_idx==0
                ref_idx=find_tetrahedron_containing_point(mesh,x_ref)
            end

        else
            make=[]
        end

        for opr in make
            V=ComplexF64[]
            I=UInt32[] #TODO: consider preallocation
            J=UInt32[]
            if opr==:M
                matrix=true #sentinel value to toggle assembley of matrix (true) or vector (false)
                for smplx in tetrahedra[mesh.domains[domain]["simplices"]]
                    CT=CooTrafo(mesh.points[:,smplx[1:4]])
                    ii,jj=create_indices(smplx)
                    vv=mass(CT)
                    append!(V,vv[:])
                    append!(I,ii[:])
                    append!(J,jj[:])
                    if output
                        ProgressMeter.next!(p)
                    end
                end
                func=(pow2,)
                arg=((:ω,),)
                txt="ω^2"
                mat="M"
            elseif opr==:K
                matrix=true
                for (smplx,c) in zip(tetrahedra[mesh.domains[domain]["simplices"]],C_tet[mesh.domains[domain]["simplices"]])
                    CT=CooTrafo(mesh.points[:,smplx[1:4]])
                    ii,jj=create_indices(smplx)
                    #println("#############")
                    #println("$CT")
                    vv=stiff(CT,c)
                    append!(V,vv[:])
                    append!(I,ii[:])
                    append!(J,jj[:])
                    if output
                        ProgressMeter.next!(p)
                    end
                end
                func=stiff_func
                arg=stiff_arg
                txt=stiff_txt
                #V=-V
                mat="K"
            elseif opr==:C
                matrix=true
                for (smplx,c) in zip(triangles[mesh.domains[domain]["simplices"]],C_tri[mesh.domains[domain]["simplices"]])
                    CT=CooTrafo(mesh.points[:,smplx[1:3]])
                    ii,jj=create_indices(smplx)
                    vv=bound(CT,c)
                    append!(V,vv[:])
                    append!(I,ii[:])
                    append!(J,jj[:])
                    if output
                        ProgressMeter.next!(p)
                    end
                end
                V.*=-1im
                func=boundary_func # func=(pow1,pow1,)
                arg=boundary_arg  # arg=((:ω,),(adm_sym,),)
                txt=boundary_txt  # txt=adm_txt
                mat="C"
            elseif opr==:Q
                matrix=true
                S=ComplexF64[]
                G=ComplexF64[]
                for smplx in tetrahedra[mesh.domains[domain]["simplices"]]
                    CT=CooTrafo(mesh.points[:,smplx[1:4]])
                    mm=volsrc(CT)
                    append!(S,mm[:])
                    append!(I,smplx[:])
                    if output
                        ProgressMeter.next!(p)
                    end
                end
                smplx=tetrahedra[ref_idx]
                CT=CooTrafo(mesh.points[:,smplx[1:4]])
                mm=gradsrc(CT,n_ref,x_ref)
                append!(G,mm[:])
                append!(J,smplx[:])
                G=-nlocal.*G
                I,J,V=outer(I,S,J,G)
                func =flame_func
                arg=flame_arg
                txt=flame_txt
                mat="Q"
            elseif opr==:m
                matrix=false
                for (smplx,c) in zip(triangles[mesh.domains[domain]["simplices"]],C_tri[mesh.domains[domain]["simplices"]])
                    CT=CooTrafo(mesh.points[:,smplx[1:3]])
                    ii=smplx[:]
                    vv=wallsrc(CT,c)
                    append!(V,vv[:])
                    append!(I,ii[:])
                    if output
                        ProgressMeter.next!(p)
                    end
                end
                V./=1im
                func=(boundary_func...,pow1,)
                arg=(boundary_arg...,(speak_sym,),)
                txt="speaker"
                mat="m"
            end

            if matrix
                #assemble discretization matrices in sparse format
                if bloch
                    for (i,j,v,f,a,t) in zip(blochify(I,J,V,naxis,nxbloch,nsector,naxis_ln,nsector_ln,N_points)...,[(),(exp_plus,),(exp_minus,),(bloch_filt,),(bloch_exp_plus,),(bloch_exp_minus,),],[(), ((b,),), ((b,),),((b,),),((b,),),((b,),),], ["", txt_plus, txt_minus, txt_filt, txt_filt_plus, txt_filt_minus,])
                        M =SparseArrays.sparse(i,j,v,dim,dim)
                        push!(L,Term(M,(func...,f...),(arg...,a...),txt*t,mat))
                    end
                else
                    M =SparseArrays.sparse(I,J,V,dim,dim)
                    push!(L,Term(M,func,arg,txt,mat))
                end
            else
                #assemble discretization vectors in sparse format
                M=SparseArrays.sparsevec(I,V,dim)
                push!(rhs,Term(M,func,arg,txt,mat))
            end

        end
    end


    if mass_weighting||bloch
        V=ComplexF64[]
        I=UInt32[] #IDEA: consider preallocation
        J=UInt32[]
        for smplx in tetrahedra
            CT=CooTrafo(mesh.points[:,smplx[1:4]])
            ii,jj=create_indices(smplx)
            vv=mass(CT)
            append!(V,vv[:])
            append!(I,ii[:])
            append!(J,jj[:])
        end
    end
    if bloch

        #IDEA: implement switch to select weighting matrix

        IJV=blochify(I,J,V,naxis,nxbloch,nsector,naxis_ln,nsector_ln,N_points,axis=false)#TODO check whether this matrix needs more blochfication... especially because tetrahedra is not periodic
        I=[IJV[1][1]..., IJV[1][2]..., IJV[1][3]...]
        J=[IJV[2][1]..., IJV[2][2]..., IJV[2][3]...]
        V=[IJV[3][1]..., IJV[3][2]..., IJV[3][3]...]
        M=SparseArrays.sparse(I,J,-V,dim,dim)

        if 0<naxis #modify axis dof for essential boundary condition when blochwave number !=0
            DI=1:naxis
            DV=ones(ComplexF64,naxis) #NOTE: make this more compact
            if order==:quad
                DI=vcat(DI,(N_points+1:naxis_ln).-nxbloch)
                DV=vcat(DV,ones(ComplexF64,naxis_ln-N_points))
            end
            for idx=1:naxis
                DV[idx]=1/M[idx,idx]
            end
            if order==:quad
                for (idx,ii) in enumerate((N_points+1:naxis_ln).-nxbloch)
                    DV[idx+naxis]=1/M[ii,ii]
                end
            end
            DM=SparseArrays.sparse(DI,DI,DV,dim,dim)
            push!(L,Term(DM,(anti_bloch_filt,),((:b,),),"(1-δ(b))","D"))
        end
    end

    if mass_weighting&&!bloch
        M=SparseArrays.sparse(I,J,-V,dim,dim)
    end
    push!(L,Term(M,(pow1,),((:λ,),),"-λ","__aux__"))

    if source #experimental return mode returning the source term
        return L,rhs
    else #classic return mode
        return L
    end
end



end #module Helmholtz
##
# module lin
# export assemble_volume_source, assemble_gradient_source, assemble_mass_operator, assemble_stiffness_operator, assemble_boundary_mass_operator
# include("FEMlin.jl")
# end
# module quad20
# export assemble_volume_source, assemble_gradient_source, assemble_mass_operator, assemble_stiffness_operator, assemble_boundary_mass_operator
# include("FEMquadlin.jl")
# end
#
# module quad
# export assemble_volume_source, assemble_gradient_source, assemble_mass_operator, assemble_stiffness_operator, assemble_boundary_mass_operator
# include("FEMquad.jl")
# end
#
# import .lin
# import .quad
# import .quad20
#
# """
#     L=discretize(mesh, dscrp, C; el_type=1, c_type=0, b=:__none__, mass_weighting=true)
#
# Discretize the Helmholtz equation using the mesh `mesh`.
#
# # Arguments
# - `mesh::Mesh`: tetrahedral mesh
# - `dscrp::Dict `: dictionary containing information on where to apply physical constraints. Such as standard wave propagation, boundary conditions, flame responses, etc.
# - `C:Array`: array defining the speed of sound. If `c_type==0` the speed of sound is constant along one tetrahedron and `length(C)==length(mesh.tetrahedra)`. If `c_type==1` the speed of sound is linearly interpolated between the vertices of the mesh and `length(C)==size(mesh.points,2)`.
# - `el_type = 1`: optional paramater to select between first (`el_type==1` the default) and second order (`el_type==2`) finite elements.
# - `c_type = 1`: optional parameter controlling whether speed of sound is constant on a tetrahedron or linearly interpolated between vertices.
# - `b::Symbol=:__none__`: optional parameter defining the Bloch wave number. If `b=:__none__` (the default) no Blochwave formalism is applied.
# - `mass_weighting=true`: optional parameter if true mass matrix is used as weighting matrix for householder, otherwise this matrix is not set.
#
# # Returns
# - `L::LinearOperatorFamily`: parametereized discretization of the specified Helmholtz equation.
# """
# function discretize(mesh::Mesh, dscrp, C; el_type=1, c_type=0, b=:__none__, mass_weighting=true)
#     if el_type==1 && c_type==0
#         #println("using linear finite elements and locally constant speed of sound")
#         assemble_volume_source, assemble_gradient_source, assemble_mass_operator, assemble_stiffness_operator, assemble_boundary_mass_operator =lin.assemble_volume_source, lin.assemble_gradient_source, lin.assemble_mass_operator, lin.assemble_stiffness_operator, lin.assemble_boundary_mass_operator
#     elseif el_type==2 && c_type==0
#         #println("using linear finite elements and locally constant speed of sound")
#         assemble_volume_source, assemble_gradient_source, assemble_mass_operator, assemble_stiffness_operator, assemble_boundary_mass_operator =quad20.assemble_volume_source, quad20.assemble_gradient_source, quad20.assemble_mass_operator, quad20.assemble_stiffness_operator, quad20.assemble_boundary_mass_operator
#     elseif el_type==2 && c_type==1
#         #println("using quadratic finite elements and locally linear speed of sound")
#         assemble_volume_source, assemble_gradient_source, assemble_mass_operator, assemble_stiffness_operator, assemble_boundary_mass_operator =quad.assemble_volume_source, quad.assemble_gradient_source, quad.assemble_mass_operator, quad.assemble_stiffness_operator, quad.assemble_boundary_mass_operator
#     else
#         println("ERROR: Elements not supported!")
#         return
#     end
#
#     L=LinearOperatorFamily(["ω","λ"],complex([0.,Inf]))
#     N_points=size(mesh.points)[2]
#     dim=deepcopy(N_points)
#     #Prepare Bloch wave analysis
#     if b!=:__none__
#         bloch=true
#         naxis=mesh.dos.naxis
#         nxbloch=mesh.dos.nxbloch
#         nsector=naxis+mesh.dos.nxsector
#         naxis_ln=mesh.dos.naxis_ln+N_points
#         nsector_ln=mesh.dos.naxis_ln+mesh.dos.nxsector_ln+N_points
#         Δϕ=2*pi/mesh.dos.DOS
#         exp_plus(z,k)=exp_az(z,Δϕ*1.0im,k) #check whether this is a closure
#         exp_minus(z,k)=exp_az(z,-Δϕ*1.0im,k)
#         bloch_filt=zeros(ComplexF64,mesh.dos.DOS)
#         bloch_filt[1]=(1.0+0.0im)/mesh.dos.DOS
#         bloch_filt=fft(bloch_filt)
#         bloch_filt=generate_Σy_exp_ikx(bloch_filt)
#         anti_bloch_filt=generate_1_gz(bloch_filt)
#         bloch_exp_plus=generate_gz_hz(bloch_filt,exp_plus)
#         bloch_exp_minus=generate_gz_hz(bloch_filt,exp_minus)
#         txt_plus="*exp(i$(b)2π/$(mesh.dos.DOS))"
#         txt_minus="*exp(-i$(b)2π/$(mesh.dos.DOS))"
#         txt_filt="*δ($b)"
#         txt_filt_plus=txt_filt*txt_plus
#         txt_filt_minus=txt_filt*txt_minus
#         #txt_filt="*1($b)"
#         #bloch_filt=pow0
#         L.params[b]=0.0+0.0im
#         dim-=mesh.dos.nxbloch
#         #N_points_dof=N_points-mesh.dos.nxbloch
#     else
#         bloch=false
#         naxis=nsector=0
#     end
#
#
#     #aggregator to create local elements ( tetrahedra, triangles)
#     if el_type==2
#         triangles=Array{UInt32,1}[] #TODO: preallocation?
#         tetrahedra=Array{UInt32,1}[]
#         tet=Array{UInt32}(undef,10)
#         tri=Array{UInt32}(undef,6)
#         for (idx,smplx) in enumerate(mesh.tetrahedra) #TODO: no enumeration
#             tet[1:4]=smplx[:]
#             tet[5]=get_line_idx(mesh,smplx[[1,2]])+N_points#find_smplx(mesh.lines,smplx[[1,2]])+N_points #TODO: type stability
#             tet[6]=get_line_idx(mesh,smplx[[1,3]])+N_points
#             tet[7]=get_line_idx(mesh,smplx[[1,4]])+N_points
#             tet[8]=get_line_idx(mesh,smplx[[2,3]])+N_points
#             tet[9]=get_line_idx(mesh,smplx[[2,4]])+N_points
#             tet[10]=get_line_idx(mesh,smplx[[3,4]])+N_points
#             push!(tetrahedra,copy(tet))
#         end
#         for (idx,smplx) in enumerate(mesh.triangles)
#             tri[1:3]=smplx[:]
#             tri[4]=get_line_idx(mesh,smplx[[1,2]])+N_points
#             tri[5]=get_line_idx(mesh,smplx[[1,3]])+N_points
#             tri[6]=get_line_idx(mesh,smplx[[2,3]])+N_points
#             push!(triangles,copy(tri))
#         end
#         dim+=size(mesh.lines)[1]
#         if bloch
#             dim-=mesh.dos.nxbloch_ln
#         end
#
#     elseif el_type==1
#             tetrahedra=mesh.tetrahedra
#             triangles=mesh.triangles
#     end
#
#     if c_type==0
#         C_tet=C
#         if mesh.tri2tet[1]==0xffffffff
#             link_triangles_to_tetrahedra!(mesh)
#         end
#         C_tri=C[mesh.tri2tet]
#     elseif c_type==1
#         C_tet=Array{UInt32,1}[] #TODO: Preallocation
#         C_tri=Array{UInt32,1}[]
#         for tet in tetrahedra
#             push!(C_tet,C[tet[1:4]])
#         end
#         for tri in triangles
#             push!(C_tri,C[tri[1:3]])
#         end
#     end
#
#     ## build matrices from domain definitions
#     for (domain,(type,data)) in dscrp
#         if type==:interior
#             make=[:M,:K]
#
#         elseif type==:admittance
#             make=[:C]
#             if length(data)==2
#                 adm_sym,adm_val=data
#
#                 adm_txt="ω*"*string(adm_sym)
#                 if adm_sym ∉ keys(L.params)
#                     L.params[adm_sym]=adm_val
#                 end
#                 boundary_func=(pow1,pow1,)
#                 boundary_arg=((:ω,),(adm_sym,),)
#                 boundary_txt=adm_txt
#             elseif length(data)==1
#                 boundary_func=(generate_z_g_z(data[1]),)
#                 boundary_arg=((:ω,),)
#                 boundary_txt="ω*Y(ω)"
#             elseif length(data)==4
#                 Ass,Bss,Css,Dss = data
#                 adm_txt="ω*C_s(iωI-A)^{-1}B"
#                 func_z_stsp = generate_z_g_z(generate_stsp_z(Ass,Bss,Css,Dss))
#                 boundary_func=(func_z_stsp,)
#                 boundary_arg=((:ω,),)
#                 boundary_txt=adm_txt
#             end
#
#         elseif type==:flame #TODO: unified interface with flameresponse and normalize n_ref
#             make=[:Q]
#             gamma,rho,nglobal,x_ref,n_ref,n_sym,tau_sym,n_val,tau_val=data
#             nlocal=(gamma-1)/rho*nglobal/compute_size!(mesh,domain)
#             if n_sym ∉ keys(L.params)
#                 L.params[n_sym]=n_val
#             end
#             if tau_sym ∉ keys(L.params)
#                 L.params[tau_sym]=tau_val
#             end
#             flame_func=(pow1,exp_delay,)
#             flame_arg=((n_sym,),(:ω,tau_sym))
#             flame_txt="$(string(n_sym))*exp(-iω$(string(tau_sym)))"
#
#         elseif type==:flameresponse
#             make=[:Q]
#             gamma,rho,nglobal,x_ref,n_ref,eps_sym,eps_val=data
#             nlocal=(gamma-1)/rho*nglobal/compute_size!(mesh,domain)
#             if eps_sym ∉ keys(L.params)
#                 L.params[eps_sym]=eps_val
#             end
#             flame_func=(pow1,)
#             flame_arg=((eps_sym,),)
#             flame_txt="$(string(eps_sym))"
#         elseif type==:fancyflame
#             make=[:Q]
#             gamma,rho,nglobal,x_ref,n_ref,n_sym,tau_sym, a_sym, n_val, tau_val, a_val=data
#             nlocal=(gamma-1)/rho*nglobal/compute_size!(mesh,domain)
#             if n_sym ∉ keys(L.params)
#                 L.params[n_sym]=n_val
#             end
#             if tau_sym ∉ keys(L.params)
#                 L.params[tau_sym]=tau_val
#             end
#             if a_sym ∉ keys(L.params)
#                 L.params[a_sym]=a_val
#             end
#
#             flame_func=(pow1,exp_ax2mxit,)
#             flame_arg=((n_sym,),(:ω,tau_sym,a_sym,))
#             flame_txt="$(string(n_sym))* exp($(string(a_sym))ω^2-iω$(string(tau_sym)))"
#         else
#             make=[]
#         end
#
#         for opr in make
#             if opr==:M
#                 I,J,V=assemble_mass_operator(mesh.points,tetrahedra[mesh.domains[domain]["simplices"]])
#                 func=(pow2,)
#                 arg=((:ω,),)
#                 txt="ω^2"
#                 mat="M"
#             elseif opr==:K
#                 I,J,V=assemble_stiffness_operator(mesh.points,tetrahedra[mesh.domains[domain]["simplices"]],C_tet[mesh.domains[domain]["simplices"]])
#                 func=()
#                 arg=()
#                 txt=""
#                 #V=-V
#                 mat="K"
#             elseif opr==:C
#                 I,J,V=assemble_boundary_mass_operator(mesh.points,triangles[mesh.domains[domain]["simplices"]],C_tri[mesh.domains[domain]["simplices"]])
#                 V*=-1im
#                 func=boundary_func # func=(pow1,pow1,)
#                 arg=boundary_arg  # arg=((:ω,),(adm_sym,),)
#                 txt=boundary_txt  # txt=adm_txt
#                 mat="C"
#             elseif opr==:Q
#                 I,S=assemble_volume_source(mesh.points,tetrahedra[mesh.domains[domain]["simplices"]])
#                 ref_idx=find_tetrahedron_containing_point(mesh,x_ref)
#                 J,G=assemble_gradient_source(mesh.points,tetrahedra[ref_idx],x_ref,n_ref)
#                 G=-nlocal*G
#                 #G=-G
#                 I,J,V=outer(I,S,J,G)
#                 func =flame_func
#                 arg=flame_arg
#                 txt=flame_txt
#                 mat="Q"
#             end
#
#             if bloch
#                 for (i,j,v,f,a,t) in zip(blochify(I,J,V,naxis,nxbloch,nsector,naxis_ln,nsector_ln,N_points)...,[(),(exp_plus,),(exp_minus,),(bloch_filt,),(bloch_exp_plus,),(bloch_exp_minus,),],[(), ((b,),), ((b,),),((b,),),((b,),),((b,),),], ["", txt_plus, txt_minus, txt_filt, txt_filt_plus, txt_filt_minus,])
#                     M =SparseArrays.sparse(i,j,v,dim,dim)
#                     push!(L,Term(M,(func...,f...),(arg...,a...),txt*t,mat))
#                 end
#
#             else
#                 M =SparseArrays.sparse(I,J,V,dim,dim)
#                 push!(L,Term(M,func,arg,txt,mat))
#             end
#
#         end
#     end
#
#
#     if mass_weighting||bloch
#         I,J,V=assemble_mass_operator(mesh.points,tetrahedra)
#     end
#     if bloch
#
#         #TODO: implement switch to select weighting matrix
#
#         IJV=blochify(I,J,V,naxis,nxbloch,nsector,naxis_ln,nsector_ln,N_points,axis=false)#TODO check whether this matrix needs more blochfication... especially because tetrahedra is not periodic
#         I=[IJV[1][1]..., IJV[1][2]..., IJV[1][3]...]
#         J=[IJV[2][1]..., IJV[2][2]..., IJV[2][3]...]
#         V=[IJV[3][1]..., IJV[3][2]..., IJV[3][3]...]
#         M=SparseArrays.sparse(I,J,-V,dim,dim)
#
#         if 0<naxis #modify axis dof for essential boundary condition when blochwave number !=0
#             DI=1:naxis
#             DV=ones(ComplexF64,naxis) #TODO: make this more compact
#             if el_type==2
#                 DI=vcat(DI,(N_points+1:naxis_ln).-nxbloch)
#                 DV=vcat(DV,ones(ComplexF64,naxis_ln-N_points))
#             end
#             for idx=1:naxis
#                 DV[idx]=1/M[idx,idx]
#             end
#             if el_type==2
#                 for (idx,ii) in enumerate((N_points+1:naxis_ln).-nxbloch)
#                     DV[idx+naxis]=1/M[ii,ii]
#                 end
#             end
#             DM=SparseArrays.sparse(DI,DI,DV,dim,dim)
#             push!(L,Term(DM,(anti_bloch_filt,),((:b,),),"(1-δ(b))","D"))
#         end
#     end
#
#     if mass_weighting&&!bloch
#         M=SparseArrays.sparse(I,J,-V,dim,dim)
#     end
#     push!(L,Term(M,(pow1,),((:λ,),),"-λ","__aux__"))
#
#     return L
# end
#end#module
