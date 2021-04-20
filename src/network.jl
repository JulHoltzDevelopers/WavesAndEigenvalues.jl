"""
    Network

Tools for creating axial (i.e., 1d) thermoacoustic network models.
The formulation is Riemann-based:
p=F exp(+iω*l/c)+G exp(-iω*l/c)
Au=A/ρc[ F exp(+iω*l/c)-G exp(-iω*l/c)]
"""
module Network
using WavesAndEigenvalues.NLEVP
export discretize
include("NLEVP_exports.jl")


"""
    out=duct(l,c,A,ρ=1.4*101325/c^2)

Create local matrix coefficients for duct element.

# Arguments
- `l`: length of the duct
- `c`: speed of sound
- `A` cross-sectional area of the duct
- `ρ` density
"""
function duct(l,c,A,ρ=1.4*101325/c^2)
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



"""
    out=terminal(R,c,A,ρ=1.4*101325/c^2;init=true)

Create local matrix coefficients for terminal element.

# Arguments
- `R`: reflection coefficient
- `c`: speed of sound
- `A`: cross-sectional area of the terminal
- `ρ`: density
- `init`: (optional) if true its the left end else it is the right one.
"""
function terminal(R,c,A,ρ=1.4*101325/c^2;init=true)
    #if hasmethod(R,(ComplexF64,Int))
    if typeof(R)<:Number
    else
        println("Error, unknown impedance type $(typeof(Z))")
    end

    if init
         M=[+R      -1;
           +1       +1;
           1*A/(ρ*c)  -1*A/(ρ*c)]

    else
        M=[ -1       -1;
            -1*A/(ρ*c)  +1*A/(ρ*c);
          -1     +R]
    end
    return [[M,(),()],]
end

"""
    out=flame(c1,c2,A,ρ=1.4*101325/c1^2)

Create local matrix coefficients for flame element with n-tau dynamics.

# Arguments
- `c1`: speed of sound on the unburnt side
- `c2`: speed of sound on the burnt side
- `A`: cross-sectional area of the terminal
- `ρ`: density
"""
function flame(c1,c2,A,ρ=1.4*101325/c1^2)
    M= [0        0;
        0        0;
        0         0;
        1         -1]

    out=duct(0,c1,A,ρ)
    push!(out,[M*(c2^2/c1^2-1)*A/(ρ*c1), (pow1,exp_delay), ((:n,),(:ω,:τ),),])
    return out
end

"""
    out=helmholtz(V,l_n,d_n,c,A,ρ=1.4*101325/c1^2)

Create local matrix coefficients for a sidewall Helmholtz damper element.

# Arguments
- `V`: Volume of the damper
- `l_n`: length of the damper neck
- `d_n`: diameter of the damper neck
- `c`: speed of sound
- `A`: cross-sectional area of the terminal
- `ρ`: density

# Notes:

The model is build on  the transfer matrix from [1]

# References

[1] F. P. Mechel, Formulas of Acoustic, Springer, 2004, p. 728
"""
function helmholtz(V,l_n,d_n,c,A,ρ=1.4*101325/c^2)
    #=====

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
            +1.0im*pow1(ω,k)*l/S_n-1.0im*c^2/V*pow(ω,k,-1))
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
    #push!(out,[M21, (admittance,), ((:ω,),),])
    push!(out,[-M21/ρ, (admittance,), ((:ω,),),])
    return out
end


"""
    sidewallimp(imp,c,A,ρ=1.4*101325/c^2)

Use function predefined  function `imp` to prescripe a frequency-dependent
side wall impedance.
"""
function sidewallimp(imp,c,A,ρ=1.4*101325/c1^2)
    M21= [0        0;
        -1       -1;
        0         0;
        0         0]
    function admittance(ω::ComplexF64,k::Int)::ComplexF64
        let Z=imp
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
    #push!(out,[M21*A*(ρ*c), (admittance,), ((:ω,),),])
    push!(out,[M21, (admittance,), ((:ω,),),])
end
"""
    out=lhr(V,l_n,d_n,c,A,ρ=1.4*101325/c^2)

Create linear Helmholtz model according to [1].

#Reference
[1] Nonlinear effects in acoustic metamaterial based on a cylindrical pipe with
    ordered Helmholtz resonators,J. Lan, Y. Li, H. Yu, B. Li, and X. Liu, Phys.
    Rev. Lett. A., 2017, [doi:10.1016/j.physleta.2017.01.036](https://doi.org/10.1016/j.physleta.2017.01.036)
"""
function lhr(V,l_n,d_n,c,A,ρ=1.4*101325/c^2)
    r_n=d_n/2
    S_n=pi*r_n^2
    B0=ρ*c^2
    η=1.5E-5 #dynamic viscosity or air in m^2/s
    R_vis=ρ*l_n/r_n*sqrt(η/2)*S_n# *ω^(1/2) #!!!! wrong in Lan it must be divided not miultiplied by r_n
    R_rad=1/4*ρ*r_n^2/c*S_n#*ω^2
    l=l_n+1.7*r_n #length correction
    Cm=V/(ρ*c^2*S_n^2)
    Mm=ρ*l*S_n
    #δ=Rm/Mm
    ω0=1/(sqrt(Cm*Mm))
    println("M: $Mm, C: $Cm, freq: $ω0")
    C=B0*S_n/(1im*ω0^2*V)/(S_n)  #last division to convert between impedance and flow impedance
    function impedance(ω::ComplexF64,k::Int)::ComplexF64
        let C=C, R_vis=R_vis, R_rad=R_rad, ω0=ω0
            return C*pow1(ω,k)-C*ω0^2*pow(ω,k,-1)-C*1im*R_vis/Mm*pow(ω,k,1/2)-C*1im*R_rad/Mm*pow2(ω,k)
        end
    end

    return sidewallimp(impedance,c,A,ρ)
end


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
This is an example for the discretization of a Rijke tube, with a Helmholtz
damper, a velocity node at the inlet and a pressure node at the outlet:
    network=[(:unode,(347.0, 0.01), ),
            (:duct,(0.25, 347.0, 0.01)),
            (:flame,(347,347*2,0.01)),
            (:duct,(0.25, 347.0*2, 0.01)),
            (:helmholtz,(0.02^3,0.01,.005,347*2,0.01)),
            (:duct,(0.25, 347.0*2, 0.01)),
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
        elseif elmt==:helmholtz2
            terms=lhr(data...)
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
