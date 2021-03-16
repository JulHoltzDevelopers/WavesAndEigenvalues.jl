"""
    Network

Tools for creating axial (i.e., 1d) thermoacoustic network models.
The formulation is Riemann based:
p=F exp(+iω*l/c)+G exp(-iω*l/c)
Au=A/ρc[ F exp(+iω*l/c)-G exp(-iω*l/c)]
"""
module Network
using WavesAndEigenvalues.NLEVP
export discretize
include("NLEVP_exports.jl")


"""
    out=duct(l,c,A,ρ)

Create local matrix coefficients for duct element.
"""
function duct(l,c,A,ρ)
    M=[-1      -1;        #p_previous-p_current=0
        -A/(ρ*c) A/(ρ*c); #A u_previous-A u_current=0
        0          0;     #p_current-p_next=0
        0          0]     #A u_current-A u_next=0

    M31=[0        0;
        0         0;
        1         0;
        0         0]

    M32= [0        0;
        0         0;
        0         1;
        0         0]
    M41= [0        0;
        0        0;
        0         0;
        1         0]
    M42= [0        0;
        0         0;
        0         0;
        0         -1]
    out=[[M,(),()],
         [M31, (generate_exp_az(1im*l/c),), ((:ω,),),],
         [M32, (generate_exp_az(-1im*l/c),), ((:ω,),),],
         [M41*A/(ρ*c), (generate_exp_az(1im*l/c),), ((:ω,),),],
         [M42*A/(ρ*c), (generate_exp_az(-1im*l/c),), ((:ω,),),],
         ]

    return out
end
duct(l,c,A)=duct(l,c,A,1.4*101325/c^2)



function terminal(R,c,A,ρ;init=true)
    #if hasmethod(R,(ComplexF64,Int))
    if typeof(R)<:Number
    else
        println("Error, unknown impedance type $(typeof(Z))")
    end

    if init
         M=[+1      -R;
           +1       +1;
           1*A/(ρ*c)  -1*A/(ρ*c)]

    else
        M=[ -1       -1;
            -1*A/(ρ*c)  +1*A/(ρ*c);
          -R     +1]
    end
    return [[M,(),()],]
end
terminal(R,c,A;init=init)=terminal(R,c,A,1.4*101325/c^2;init)

function flame(c1,c2,A,ρ)
    M= [0        0;
        0        0;
        0         0;
        1         -1]

    out=duct(0,c1,A,ρ)
    push!(out,[M*(c2^2/c1^2-1)*A/(ρ*c1), (pow1,exp_delay), ((:n,),(:ω,:τ),),])
    return out
end
flame(c1,c2,A)=flame(c1,c2,A,1.4*101325/c1^2)

function helmholtz(V,l_n,d_n,c,A,ρ)
    #=====
    Transfer matrix from Mechel (Formulas of Acoustic p. 728):
    p_u=p_d
    u_u=1/Z p_d+u_d <==> A u_u=A/Z p_d+A u_d>

    TODO:
    <==> Z*A*u_u=A*p_d+Z*Au_d (roe easier representation)

    where

    Z(ω)=ρ[ω^2/(π*c)*[2-r_n/r_u]+0.425*M*c*/S_n+1im[ω*l/S_n-c^2/(ω*V)]]

    where
    ρ: density
    c: speed of sound
    r_n: neck radius of the Helmholtz damper
    r_u: radius of the main duct
    M: Mach number
    S_n = π*r_n2 : crossectional area of the neck
    V = : volume of the helmholtz damper
    and
    l=l_n+t_w+0.85 r_n*(2-r_n/r_u)

    with tw: wall thickness?

    ## adaption to Riemann invariants

    p_u-A_+*exp(im*ω*l/c)+A_-*exp(-im*ω*l/c)=0
    u_u-A_+/ρc*exp(im*ω*l/c)+A_+/ρc*exp(-im*ω*l/c)=0
    The last equation is to be multiplied by A.

    ## notes
    this element does not resolve what is between the usptream and the downstream
    site. It directly describes the jump!
    =====#

    r_n=d_n/2
    r_u=sqrt(A/pi)
    S_n=pi*r_n^2
    t_w=0
    l=l_n+t_w+0.85*r_n*(2-r_n/r_u)
    M=0
    M21= [0        0;
        -1       -1;
        0         0;
        0         0]
    #M22= [0        0;
    #    0         -1;
    #    0         0;
    #    0         0]

    function impedance(ω::ComplexF64,k::Int)::ComplexF64
        let ρ=ρ, r_n=r_n, r_u=r_u, M=M, c=c, l=l, V=V ,S_n=S_n
            return ρ*(pow2(ω,k)/(π*c)*(2-r_n/r_u)+pow0(ω,k)*0.425*M*c/S_n
            +1.0im*pow1(ω,k)*l/S_n-1.0im-c^2/V*pow(ω,k,-1))
        end
    end
    function admittance(ω::ComplexF64,k::Int)::ComplexF64
        let Z=impedance
            if k==0
                return 1/Z(ω,0)
            elseif k==1
                return -Z(ω,1)/Z(ω,0)^2 #TODO: implement arbitray order chain rule
            else
                return NaN+NaN*im
            end
        end
    end
    function admittance(ω::Symbol)::String
        return "1/Z($ω)"
    end


    out=duct(0,c,A,ρ)
    push!(out,[M21*A, (admittance,), ((:ω,),),])
    #push!(out,[M22*A, (admittance,), ((:ω,),),])
end
helmholtz(V,l_n,d_n,c,A)=helmholtz(V,l_n,d_n,c,A,1.4*101325/c^2)


"""
    L=discretize(network)

Discretize network model.

# Arguments
- `network`

# Returns
- `L::LinearOperatorFamily`

# The network is specified of a list of network elements.
elements together with there options. Available elements are.
- `:duct`: duct element with options `(l,c,A,ρ)`
- `:flame`: n-tau-flame element with options `(c1,c2,A,ρ)`
- `:helmholtz`:Helmholtz damper element with options `(V,l_n,d_n,c,A,ρ)`
- `:unode`: sound hard boundary, i.e. a velocity node with options `(c,A,ρ)`
- `:pnode`: sound soft boundary, i.e., a pressure node with options `(c,A,ρ)`
- `:anechoic`: anechoic boundary condition with option `(c,A,ρ)`

For all elements the density `ρ` is an optional parameter which is computed
as `ρ=1.4*101325/c^2` if unspecified (air at atmospheric pressure).

# Example
This is an example for the discretization of a Rijke tube.
    network=[(:unode,(347.0, 0.01), ),
            (:duct,(0.25, 347.0, 0.01)),
            (:flame,(347,347*2,0.01)),
            (:duct,(0.25, 347.0*2, 0.01)),
            (:helmholtz,(0.02^3,0.01,.005,347*2,0.01)),
            (:pnode,(347.0*2, 0.01), )
            ]

    L=discretize(network)
    L.params[n]=1
    L.params[τ]=0.001
"""
function discretize(network)
    i,j=1,1
    L=LinearOperatorFamily((:ω,))
    N=length(network)
    A=zeros(2*N,2*N)
    for (idx,(elmt, data)) in enumerate(network)
        if elmt==:duct
            terms=duct(data...)
        elseif elmt==:unode
            R=1
            if idx==1
                init=true
            elseif idx==length(network)
                init=false
            else
                println("Error, terminal element at intermediate position $idx!")
            end
            #terms=velocity_node(data...)
            terms=terminal(R,data...,init=init)
        elseif elmt==:pnode
            R=-1
            if idx==1
                init=true
            elseif idx==length(network)
                init=false
            else
                println("Error, terminal element at intermediate position $idx!")
            end
            terms=terminal(R,data...,init=init)
        elseif elmt==:anechoic
            R=0
            if idx==1
                init=true
            elseif idx==length(network)
                init=false
            else
                println("Error, terminal element at intermediate position $idx!")
            end
            terms=terminal(R,data...,init=init)
        elseif elmt==:pnodeend
            terms=pressure_node_end(data...)
        elseif elmt==:vnodestart
            terms=velocity_node_start(data...)
        elseif elmt==:flame
            terms=flame(data...)
        elseif elmt==:helmholtz
            terms=helmholtz(data...)
        else
            println("Warning: unknown element $elmt")
            continue
        end
        I,J=size(terms[1][1])

        for (coeff,func,arg) in terms
            M=deepcopy(A)
            M[i:i+I-1,j:j+J-1]=coeff
            push!(L,Term(M,func,arg,"M"))
        end
        i+=I-2
        j+=2
    end
    return L
end
end #module network

## legacy code

# function velocity_node(c,A,ρ)
#     M=[1*A/(ρ*c) -1*A/(ρ*c);
#        +1       +1;
#        1*A/(ρ*c)  -1*A/(ρ*c)]
#     return [[M,(),()],]
# end
# velocity_node(c,A)=velocity_node(c,A,1.4*101325/c^2)
#
#
# function pressure_node(c,A,ρ)
#     M=[ -1       -1;
#         -1*A/(ρ*c)  +1*A/(ρ*c);
#        +1       +1]
#
#     return [[M,(),()],]
# end
# pressure_node(c,A) = pressure_node(c,A,1.4*101325/c^2)


# function velocity_node_start(l,c,A,ρ)
#     M=[A/(ρ*c) -A/(ρ*c);
#         0          0;
#         0          0]
#
#     M21=[0         0;
#         1         0;
#         0         0]
#
#     M22= [0         0;
#         0         1;
#         0         0]
#     M31= [0         0;
#         0         0;
#         1         0]
#     M32= [0         0;
#         0         0;
#         0         -1]
#     out=[[M,(),()],
#          [M21, (generate_exp_az(1im*l/c),), ((:ω,),),],
#          [M22, (generate_exp_az(-1im*l/c),), ((:ω,),),],
#          [M31*A/(ρ*c), (generate_exp_az(1im*l/c),), ((:ω,),),],
#          [M32*A/(ρ*c), (generate_exp_az(-1im*l/c),), ((:ω,),),],
#          ]
#
#     return out
# end
# velocity_node_start(l,c,A) = velocity_node_start(l,c,A,1.4*101325/c^2)
#
#
# function pressure_node_end(l,c,A,ρ)
#     M=[-1      -1;
#         -A/(ρ*c) A/(ρ*c);
#         0          0]
#
#     M31=[0        0;
#         0         0;
#         1         0]
#
#     M32= [0        0;
#         0         0;
#         0         1]
#
#     out=[[M,(),()],
#          [M31, (generate_exp_az(1im*l/c),), ((:ω,),),],
#          [M32, (generate_exp_az(-1im*l/c),), ((:ω,),),],
#          ]
#     return out
# end
# pressure_node_end(l,c,A) = pressure_node_end(l,c,A,1.4*101325/c^2)
