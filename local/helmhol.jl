include("./src/Bloch.jl")
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


function helmhol(mesh::Mesh, dscrp, C; order=:1, b=:__none__, mass_weighting=true, source=false)
    triangles,tetrahedra,dim=aggregate_elements(mesh,order)
    N_points=size(mesh.points,2)


    if length(c)==length(mesh.tetrahedra)
        C_tet=C
        if mesh.tri2tet[1]==0xffffffff
            link_triangles_to_tetrahedra!(mesh)
        end
        C_tri=C[mesh.tri2tet]
    elseif length(c)==size(mesh.points,2)
        C_tet=Array{UInt32,1}[] #TODO: Preallocation
        C_tri=Array{UInt32,1}[]
        for tet in tetrahedra
            push!(C_tet,C[tet[1:4]])
        end
        for tri in triangles
            push!(C_tri,C[tri[1:3]])
        end
    end

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
        if order==:1
            dim-=mesh.dos.nxbloch
        elseif order==:2
            dim-=mesh.dos.nxbloch+mesh.dos.nxbloch_ln
        elseif order==:h
            dim-=4*mesh.dos.nxbloch
        end
        #N_points_dof=N_points-mesh.dos.nxbloch
    else
        bloch=false
        naxis=nsector=0
    end
    #wrapper for FEM constructors TODO: use multipledispatch and move to FEM.jl
    function stiff(J,c)
        if order==:1
            if length(c)==1
                return -c^2*s43nv1nu1(J)
            elseif length(c)==4
                return -s43nv1nu1cc1(J,c)
            end
        elseif order==:2
            if length(c)==1
                return -c^2*s43nv2nu2(J)
            elseif length(c)==4
                return -s43nv2nu2cc1(J,c)
            end
        elseif order==:h
            if length(c)==1
                return -c^2*s43nvhnuh(J)
            elseif length(c)==4
                return -s43nvhnuhcc1(J,c)
            end
        end
    end
    function mass(J)
        if order==:1
            return s43v1u1(J)
        elseif order==:2
            return s43v2u2(J)
        elseif order==:h
            return s43vhuh(J)#massh(J)
        end
    end

    function bound(J,c)
        if order==:1
            if length(c)==1
                return c*s33v1u1(J)
            elseif length(c)==4
                return s33v1u1c1(J,c)
            end
        elseif order==:2
            if length(c)==1
                return c*s33v2u2(J)
            elseif length(c)==4
                return s33v2u2c1(J,c)
            end
        elseif order==:h
            if length(c)==1
                return c*s33vhuh(J)
            elseif length(c)==4
                return s33vhuhc1(J,c)
            end
        end
    end

    function volsrc(J)
        if order==:1
            return s43v1(J)
        elseif order==:2
            return s43v2(J)
        elseif order==:h
            return s43vh(J)
        end
    end

    function gradsrc(J,n,x)
        if order==:1
            return s43nv1rx(J,n,x)
        elseif order==:2
            return s43nv2rx(J,n,x)
        elseif order==:h
            return s43nvhrx(J,n,x)
        end
    end

    function wallsrc(J)
        if order==:1
            return s33v1(J)
        elseif order==:2
            return s33v2(J)
        elseif order==:h
            return s33vh(J)
    end


    #initialize linear operator family...
    L=LinearOperatorFamily(["ω","λ"],complex([0.,Inf]))
    #...and source vector
    rhs=LinearOperatorFamily(["ω"],complex([0.,]))


    ## build discretization matrices and vectors from domain definitions
    for (domain,(type,data)) in dscrp
        if type==:interior
            make=[:M,:K]
            matrix=true #sentinel value to toggle assembley of matrix (true) or vector (false)

        elseif type==:admittance
            make=[:C]
            if length(data)==2
                adm_sym,adm_val=data

                adm_txt="ω*"*string(adm_sym)
                if adm_sym ∉ keys(L.params)
                    L.params[adm_sym]=adm_val
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

            matrix=true

        elseif type==:flame #TODO: unified interface with flameresponse and normalize n_ref
            make=[:Q]
            gamma,rho,nglobal,x_ref,n_ref,n_sym,tau_sym,n_val,tau_val=data
            nlocal=(gamma-1)/rho*nglobal/compute_size!(mesh,domain)
            if n_sym ∉ keys(L.params)
                L.params[n_sym]=n_val
            end
            if tau_sym ∉ keys(L.params)
                L.params[tau_sym]=tau_val
            end
            flame_func=(pow1,exp_delay,)
            flame_arg=((n_sym,),(:ω,tau_sym))
            flame_txt="$(string(n_sym))*exp(-iω$(string(tau_sym)))"

            matrix=true

        elseif type==:flameresponse
            make=[:Q]
            gamma,rho,nglobal,x_ref,n_ref,eps_sym,eps_val=data
            nlocal=(gamma-1)/rho*nglobal/compute_size!(mesh,domain)
            if eps_sym ∉ keys(L.params)
                L.params[eps_sym]=eps_val
            end
            flame_func=(pow1,)
            flame_arg=((eps_sym,),)
            flame_txt="$(string(eps_sym))"

            matrix=true

        elseif type==:fancyflame
            make=[:Q]
            gamma,rho,nglobal,x_ref,n_ref,n_sym,tau_sym, a_sym, n_val, tau_val, a_val=data
            nlocal=(gamma-1)/rho*nglobal/compute_size!(mesh,domain)
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
            matrix=true

        elseif type==:loudspeaker
            make=[:m]
            loud_sym,loud_val=data
            matrix=false

        else
            make=[]
        end

        for opr in make
            V=ComplexF64[]
            I=UInt32[] #TODO: consider preallocation
            J=UInt32[]
            if opr==:M
                for smplx in tetrahedra[mesh.domains[domain]["simplices"]]
                    CT=CooTrafo(mesh.points[:,smplx[1:4]])
                    ii,jj=create_indices(smplx)
                    vv=mass(CT)
                    append!(V,vv[:])
                    append!(I,ii[:])
                    append!(J,jj[:])
                end
                func=(pow2,)
                arg=((:ω,),)
                txt="ω^2"
                mat="M"
            elseif opr==:K
                for (smplx,c) in zip(tetrahedra[mesh.domains[domain]["simplices"]],C_tet[mesh.domains[domain]["simplices"]])
                    CT=CooTrafo(mesh.points[:,smplx[1:4]])
                    ii,jj=create_indices(smplx)
                    vv=stiff(CT,c) #TODO: move c^2 to stiff function
                    append!(V,vv[:])
                    append!(I,ii[:])
                    append!(J,jj[:])
                end
                func=()
                arg=()
                txt=""
                #V=-V
                mat="K"
            elseif opr==:C
                for (smplx,c) in zip(triangles[mesh.domains[domain]["simplices"]],C_tri[mesh.domains[domain]["simplices"]])
                    CT=CooTrafo(mesh.points[:,smplx[1:3]])
                    ii,jj=create_indices(smplx)
                    vv=bound(CT,c) #TODO: move c to stiff function
                    append!(V,vv[:])
                    append!(I,ii[:])
                    append!(J,jj[:])
                end
                V.*=-1im
                func=boundary_func # func=(pow1,pow1,)
                arg=boundary_arg  # arg=((:ω,),(adm_sym,),)
                txt=boundary_txt  # txt=adm_txt
                mat="C"
            elseif opr==:Q
                S=ComplexF64[]
                G=ComplexF64[]
                for smplx in tetrahedra[mesh.domains[domain]["simplices"]]
                    CT=CooTrafo(mesh.points[:,smplx[1:4]])
                    mm=volsrc(CT)
                    append!(S,mm[:])
                    append!(I,smplx[:])
                end
                ref_idx=find_tetrahedron_containing_point(mesh,x_ref)
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
                for (smplx,c) in zip(triangles[mesh.domains[domain]["simplices"]],C_tri[mesh.domains[domain]["simplices"]])
                    CT=CooTrafo(mesh.points[:,smplx[1:3]])
                    ii=smplx[:]
                    vv=wallsrc(CT,c)
                    append!(V,vv[:])
                    append!(I,ii[:])
                end
                func=(pow1,pow1,)
                arg=((:ω,),(loud_sym,),)
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
        I=UInt32[] #TODO: consider preallocation
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

        #TODO: implement switch to select weighting matrix

        IJV=blochify(I,J,V,naxis,nxbloch,nsector,naxis_ln,nsector_ln,N_points,axis=false)#TODO check whether this matrix needs more blochfication... especially because tetrahedra is not periodic
        I=[IJV[1][1]..., IJV[1][2]..., IJV[1][3]...]
        J=[IJV[2][1]..., IJV[2][2]..., IJV[2][3]...]
        V=[IJV[3][1]..., IJV[3][2]..., IJV[3][3]...]
        M=SparseArrays.sparse(I,J,-V,dim,dim)

        if 0<naxis #modify axis dof for essential boundary condition when blochwave number !=0
            DI=1:naxis
            DV=ones(ComplexF64,naxis) #TODO: make this more compact
            if el_type==2
                DI=vcat(DI,(N_points+1:naxis_ln).-nxbloch)
                DV=vcat(DV,ones(ComplexF64,naxis_ln-N_points))
            end
            for idx=1:naxis
                DV[idx]=1/M[idx,idx]
            end
            if el_type==2
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
##
mesh=Mesh("./examples/tutorials/rect_tube.msh",scale=0.001)
dscrp=Dict() #initialize model discreptor
dscrp["Interior"]=(:interior, ()) #define resonant cavity
dscrp["Outlet"]=(:admittance, (:Y,1*1E15)) #specify outlet BC
γ=1.4 #ratio of specific heats
ρ=1.225 #density at the reference location upstream to the flame in kg/m^3
Tu=300.0    #K unburnt gas temperature
Tb=1200.0    #K burnt gas temperature
P0=101325.0 # ambient pressure in Pa
A=0.025^2 # cross sectional area of the tube
Q02U0=P0*(Tb/Tu-1)*A*γ/(γ-1) #the ratio of mean heat release to mean velocity Q02U0
x_ref=[0.0; 0.0; -0.00101] #reference point
n_ref=[0.0; 0.0; 1.00] #directional unit vector of reference velocity
n=1*0.01 #interaction index
τ=0.001 #time delay
dscrp["Flame"]=(:flame,(γ,ρ,Q02U0,x_ref,n_ref,:n,:τ,n,τ)) #flame dynamics
R=287.05 # J/(kg*K) specific gas constant (air)
speedofsound(x,y,z) = z<0. ? sqrt(γ*R*Tu) : sqrt(γ*R*Tb)
c=generate_field(mesh,speedofsound)
L=helmhol(mesh,dscrp,c,order=:h)
L_val=WavesAndEigenvalues.Helmholtz.discretize(mesh,dscrp,c)


##
sol,nn,flag=householder(L,340*2*pi,maxiter=5,output=true)
##
data=Dict()
data["bla"]=real.(sol.v)
vtk_write("tube",mesh,data)
##
c=ones(length(mesh.tetrahedra))*347
c=ones(size(mesh.points,2))*347

dscrp=Dict()
dscrp["Interior"]=(:interior,())
H=WavesAndEigenvalues.Helmholtz.discretize(mesh,dscrp,c,el_type=2,c_type=1)
sol,nn,flag=householder(H,330*2*pi,maxiter=2,output=true)
###

X=[0 1 0 0;
   1 0 0 1;
   0 0 1 1.]
X=rand(3,4)*100
J=CooTrafo(X)
m1=s43vhuh(J)
m2=massh(J)
##
