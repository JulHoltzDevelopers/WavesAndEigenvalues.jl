#This file will be turned in more flexible algebra routines

##header
import SpecialFunctions

## Polynomial type
"""
Polynomial type
"""
struct Polynomial #TODO: make type parametric use eltype(f)
 coeffs
end

#make Polynomial callable with Horner-scheme evaluation
function (p::Polynomial)(z,k::Int=0)
    p=derive(p,k)
    f=0.0+0.0im #typing
    for i = length(p.coeffs):-1:1
       f*=z
       f+=p.coeffs[i]
    end
    return f
end

import Base.show
import Base.string

function string(p::Polynomial)
    N=length(p.coeffs)
    if N==0
        txt="0"
    else
        txt="$(p.coeffs[1])"
        for n=2:N
            coeff=string(p.coeffs[n])
            if occursin("+",coeff) || occursin("-",coeff)
                txt*="+($coeff)*z^$(n-1)"
            else
                txt*="+$coeff*z^$(n-1)"
            end
        end
    end
    return txt
end

function show(io::IO,p::Polynomial)
    print(io,string(p))
end



import Base.+
import Base.-
import Base.*
import Base.^

function +(a::Polynomial,b::Polynomial)
    if length(a.coeffs)>length(b.coeffs)
        a,b=b,a
    end
    c=copy(b.coeffs)
    for (idx,val) in enumerate(a.coeffs)
        c[idx]+=val
    end
    return Polynomial(c)
end

function -(a::Polynomial,b::Polynomial)
    return +(a,-b)
end

function *(a::Polynomial,b::Polynomial)
    I=length(a.coeffs)
    J=length(b.coeffs)
    K=I+J-1
    c=zeros(ComplexF64,K) #TODO: sanity check for type of b
    idx=1
    for k =1:K
        for i =1:k
            #TODO: figure out when to break the loop
            j=k-i+1
            if i>I ||j>J
                continue
            end
            c[k]+=a.coeffs[i]*b.coeffs[j]
        end
    end
    return Polynomial(c)
end

#scalar multiplication
function *(a::Polynomial,b::Number)
    return a*Polynomial([b])
end
function *(a::Number,b::Polynomial)
    return Polynomial([a])*b
end

function ^(p::Polynomial, k::Int)
    b=Polynomial(copy(p.coeffs))
    for i=1:k-1
        b*=p
    end
    return b
end

function prodrange(start,stop)
    return SpecialFunctions.factorial(stop)/SpecialFunctions.factorial(start)
end

"""
    b=derive(a::Polynomial,k::Int=1)

Compute `k`th derivative of Polynomial `a`.
"""
function derive(a::Polynomial,k::Int=1)
    N=length(a.coeffs)
    if k>=N
        return Polynomial([]) #TODO: type
    end
    d=zeros(typeof(a.coeffs[1]),N-k)
    for i=1:N-k
        d[i]=a.coeffs[i+k]*prodrange(i-1,i+k-1)
    end
    return Polynomial(d)
end

"""
    g=shift(f,Δz)

Shift Polynomial with respect to its argument. (Taylor shift)
"""
function shift(f::Polynomial,Δz)
    g=Array{eltype(f)}(undef,length(f.coeffs))
    for n=0:length(f.coeffs)-1
        g[n+1]=f(Δz)/SpecialFunctions.factorial(n*1.0)
        f=Polynomial(f.coeffs[2:end])
        f.coeffs.*=1:length(f.coeffs)
    end
    return Polynomial(g)
end

function scale(p::Polynomial,a::Number)
        p=p.coeffs
        s=1
        for idx =1: length(p)
            p[idx]*=s
            s*=a
        end
        return Polynomial(p)
end

#TODO: write macro nto create !-versions

## tests
a=Polynomial([0.0,1.0])
b=a*a
c=derive(a)
d=a+b+c
e=d^3
f=shift(e,2)
e(2)-f(0)
## rational Polynomial
struct rational_Polynomial
    num::Polynomial
    den::Polynomial
end

function derive(a::rational_Polynomial)
    return rational_Polynomial(derive(a.num)*a.den-a.num*derive(a.den),a.den^2)
end

"""
    derive(a::Any, k:Int)

Take the `k`th derivative of `a` by calling the derive function multiple times.

# Notes
This is a fallback implementation.
"""
function derive(a::Any,k::Int)
    for i=1:k
        a=derive(a)
    end
    return a
end
##
function ndd(vararg...)
    #parse input TODO: sanity check that z is unique
    #get number of points
    #(confluent points are counted multiple times)
    N=0
    for i = 1:length(vararg)÷2
        N+=length(vararg[2*i])
    end
    #initialize lists
    z=Array{ComplexF64}(undef,N) #TODO: proper type rules
    coeffs=Array{ComplexF64}(undef,N)
    cnfl=Array{Int64}(undef,N) # confluence counter array
    n=0 #total counter
    for i = 1:length(vararg)÷2
        z[1+n:n+length(vararg[2*i])] .= vararg[2*i-1]
        coeffs[1+n:n+length(vararg[2*i])]  = vararg[2*i]
        cnfl[1+n:n+length(vararg[2*i])] = 0:length(vararg[2*i])-1
        n+=length(vararg[2*i])
    end



    #pyramid scheme for computation of divided-difference scheme
    println("###########")
    println(cnfl)
    Z=copy(z)
    for i=1:N-1
        println(coeffs)
        for j=N:-1:(i+1)
            Δz=z[j]-z[j-i]
            if Δz!=0 # this is a non-confluent point
               coeffs[j]=(coeffs[j]-coeffs[j-1-cnfl[j-1]])/Δz
            end
            if cnfl[j-1]!=0
                cnfl[j-1]-=1
            end
        end
    end
    println(coeffs)


    # construct newton polynomial from divided difference coefficients
    basis=Polynomial([1])
    p=Polynomial([0])
    for i=1:N
        p+=basis*coeffs[i]
        basis*=Polynomial([-z[i],1])
    end

    return p
end


##
m0=Polynomial([1])
m1=Polynomial([-1,1])
m2=Polynomial([-2,1])
m3=Polynomial([-3,1])

p=3*m0+(-1)*m1+2.5*m1*m2
##
p=ndd(1,[3,],2,[2,],3,[6,])
p.([1,2,3])
##
p=ndd(2,[6,3,5])
p(2,3)
##
p=ndd(1,[1,2],2,[4])
p.([1,2],0)
##
vals=(1,[1,2,3],2,[4,5,0,7,8,9],3,[60,1])
p=ndd(vals...)
##
println("#########Here##########")
for i=1:length(vals)÷2
    z=vals[2*i-1]
    for (k,coeff) in enumerate(vals[2*i])
        k=k-1
        test=p(z,k)/factorial(k)-coeff==0
        if !test
            println("z:$z k:$k coeff:$coeff")
        end
    end
end

## TODO next: Find reference for this algorithm.
function newton_divided_difffunction(vararg...)
    #hässliches steigungsschema...
    coeffs=Dict()
    z=Array{ComplexF64}(undef,length(vararg)÷2)
    comb=[]
    #initialize divided difference scheme
    #use taylor coefficient at confluent points
    for (idx,n) in enumerate(1:2:length(vararg))
        z[idx]=vararg[n]
        taylor_coeffs=vararg[n+1]
        loc_tab=repeat([idx],length(taylor_coeffs))
        push!(comb,loc_tab...)
        for (i,j) in enumerate(1:length(loc_tab)) #TODO: enumeration and i are not required
            coeffs[loc_tab[1:j]]=taylor_coeffs[j]
        end
    end

    #fill missing values in the newton divided difference table
    for n =1:length(comb)
        for j =1:length(comb)-n
            idx=comb[j:j+n]
            if idx ∉ keys(coeffs)
                a=coeffs[idx[1:end-1]]
                b=coeffs[idx[2:end]]
                z_a=z[idx[1]]
                z_b=z[idx[end]]
                coeffs[idx]=(b-a)/(z_b-z_a)
            end
        end
    end

    #construct newton_Polynomial
    basis_function=Polynomial([1])
    idx=comb[1:1]


    c=Polynomial(coeffs[idx])
    p=basis_function*c
    for n=2:length(comb)
        idx=comb[1:n]
        c=Polynomial(coeffs[idx])
        basis_function=basis_function*Polynomial(-z[idx[end-1]],1]
        p=p+ᴾc*ᴾbasis_function
    end
    basis_function=basis_function*ᴾ[-z[comb[end]],1]


    #Kronecker's Algorithm to find the interpolating Pade approximant
    #initialize Kronecker's algorithm
    # (see Chapter 7.1 (page 341) in G. A. Baker and P. Graves-Morris.
    # Padé Approximants, 2nd Edition)
    p0=basis_function
    q0=[]
    q=[1]
    pade_table=[]
    pade_table=push!(pade_table,(p,q))

    while length(p)>1
        K=2 #system size

        A=zeros(ComplexF64,K,K)
        y=zeros(ComplexF64,K)
        for k=1:K
            y[k]=p0[end-k+1]
            for i=1:k
                    A[k,i]=p[end-k+i]
            end
        end
        x=A\y
        #TODO: sanity check!
        reverse!(x)
        p,p0=x*ᴾp-ᴾp0,p
        q,q0=x*ᴾq-ᴾq0,q
        p=p[1:end-K]
        push!(pade_table,(p,q))
    end

    return pade_table

end


## helper functions to directly use arrays
function *ᴾ(a,b)
    I=length(a)
    J=length(b)
    K=I+J-1
    c=zeros(ComplexF64,K)
    idx=1
    for k =1:K
        for i =1:k
            #TODO: figure out when to break the loop
            j=k-i+1
            if i>I ||j>J
                continue
            end
            c[k]+=a[i]*b[j]
        end
    end
    return c
end


function +ᴾ(a,b)
    if length(a)>length(b)
        a,b=b,a
    end
    for (idx,val) in enumerate(a)
        b[idx]+=val
    end
    return b
end

function -ᴾ(a,b)
    return +ᴾ(a,-b)
end

function derive(a)
    N=length(a)
    d=zeros(ComplexF64,N-1)
    for i=1:N-1
        d[i]=a[i+1]*i
    end
    return d
end
##
function compute_newton_polynomial(vararg...)
    #hässliches steigungsschema...
    coeffs=Dict()
    z=Array{ComplexF64}(undef,length(vararg)÷2)
    comb=[]
    #initialize divided difference scheme
    #use taylor coefficient at confluent points
    for (idx,n) in enumerate(1:2:length(vararg))
        z[idx]=vararg[n]
        taylor_coeffs=vararg[n+1]
        loc_tab=repeat([idx],length(taylor_coeffs))
        push!(comb,loc_tab...)
        for (i,j) in enumerate(1:length(loc_tab))
            coeffs[loc_tab[1:j]]=taylor_coeffs[j]
        end
    end

    #fill missing values in the newton divided difference table
    for n =1:length(comb)
        for j =1:length(comb)-n
            idx=comb[j:j+n]
            if idx ∉ keys(coeffs)
                a=coeffs[idx[1:end-1]]
                b=coeffs[idx[2:end]]
                z_a=z[idx[1]]
                z_b=z[idx[end]]
                coeffs[idx]=(b-a)/(z_b-z_a)
            end
        end
    end

    #construct newton_polynomial
    basis_function=[1]
    idx=comb[1:1]


    c=[coeffs[idx]]
    p=basis_function*ᴾc
    for n=2:length(comb)
        idx=comb[1:n]
        c=[coeffs[idx]]
        basis_function=basis_function*ᴾ[-z[idx[end-1]],1]
        p=p+ᴾc*ᴾbasis_function
    end
    basis_function=basis_function*ᴾ[-z[comb[end]],1]


    #Kronecker's Algorithm to find the interpolating pade approximant
    #initialize Kronecker's algorithm
    p0=basis_function
    q0=[]
    q=[1]
    pade_table=[]
    pade_table=push!(pade_table,(p,q))

    while length(p)>1
        K=2 #system size

        A=zeros(ComplexF64,K,K)
        y=zeros(ComplexF64,K)
        for k=1:K
            y[k]=p0[end-k+1]
            for i=1:k
                    A[k,i]=p[end-k+i]
            end
        end
        x=A\y
        #TODO: sanity check!
        reverse!(x)
        p,p0=x*ᴾp-ᴾp0,p
        q,q0=x*ᴾq-ᴾq0,q
        println("#####")
        println(p[end-K+1:end])
        p=p[1:end-K]
        push!(pade_table,(p,q))
    end

    return pade_table

end

##
polyval(p,z)= sum(p.*(z.^(0:length(p)-1)))

function taylor_shift(f,Δz)
    g=Array{eltype(f)}(undef,length(f))
    for n=0:length(f)-1
        g[n+1]=polyval(f,Δz)/SpecialFunctions.factorial(n*1.0)
        f=f[2:end]
        f.*=1:length(f)
    end
    return g
end

#Test
f=[-42.168, 2im, π, 1, 0.05]
z= π*1im-5.923779
g=taylor_shift(f,z)
polyval(f,z)-polyval(g,0)

#TODO implement Kronecker's algorithm
function multi_point_pade(L,M,vararg...;Z0=0)
    #TODO sanity checks length vararg should be divisibl #and the sum of the entries should amound to L+M

    #Build linear system



end



function newton_polynomial(vararg...)
    N=0
    for idx =1:2:length(vararg)
        N+=length(vararg[idx+1])
    end
    A=zeros(ComplexF64,N,N)
    y=zeros(ComplexF64,N)
    eqidx=1
    for idx =1:2:length(vararg)
        coeffs=ones(ComplexF64,N)
        z0=vararg[idx]
        z=[z0^i for i in 0:N-1]
        k=1
        fac=1
        for val in vararg[idx+1]
            y[eqidx]=val*fac
            A[eqidx,:]=coeffs.*z
            coeffs[k]=0
            z[k]=0
            for (i,j) in enumerate(k+1:N)
                coeffs[j]*=i
                z[j]=z0^(i-1)
            end
            fac*=k
            k+=1
            eqidx+=1
        end

    end

    x=A\y
    return x

end


function newton_divided_difffunction(vararg...)
    #hässliches steigungsschema...
    coeffs=Dict()
    z=Array{ComplexF64}(undef,length(vararg)÷2)
    comb=[]
    #initialize divided difference scheme
    #use taylor coefficient at confluent points
    for (idx,n) in enumerate(1:2:length(vararg))
        z[idx]=vararg[n]
        taylor_coeffs=vararg[n+1]
        loc_tab=repeat([idx],length(taylor_coeffs))
        push!(comb,loc_tab...)
        for (i,j) in enumerate(1:length(loc_tab))
            coeffs[loc_tab[1:j]]=taylor_coeffs[j]
        end
    end

    #fill missing values in the newton divided diffrence table
    for n =1:length(comb)
        for j =1:length(comb)-n
            idx=comb[j:j+n]
            if idx ∉ keys(coeffs)
                a=coeffs[idx[1:end-1]]
                b=coeffs[idx[2:end]]
                z_a=z[idx[1]]
                z_b=z[idx[end]]
                coeffs[idx]=(b-a)/(z_b-z_a)
            end
        end
    end

    #construct newton_polynomial
    basis_function=[1]
    idx=comb[1:1]


    c=[coeffs[idx]]
    p=basis_function*ᴾc
    for n=2:length(comb)
        idx=comb[1:n]
        c=[coeffs[idx]]
        basis_function=basis_function*ᴾ[-z[idx[end-1]],1]
        p=p+ᴾc*ᴾbasis_function
    end
    basis_function=basis_function*ᴾ[-z[comb[end]],1]


    #Kronecker's Algorithm to find the interpolating pade approximant
    #initialize Kronecker's algorithm
    p0=basis_function
    q0=[]
    q=[1]
    pade_table=[]
    pade_table=push!(pade_table,(p,q))

    while length(p)>1
        K=2 #system size

        A=zeros(ComplexF64,K,K)
        y=zeros(ComplexF64,K)
        for k=1:K
            y[k]=p0[end-k+1]
            for i=1:k
                    A[k,i]=p[end-k+i]
            end
        end
        x=A\y
        #TODO: sanity check!
        reverse!(x)
        p,p0=x*ᴾp-ᴾp0,p
        q,q0=x*ᴾq-ᴾq0,q
        println("#####")
        println(p[end-K+1:end])
        p=p[1:end-K]
        push!(pade_table,(p,q))
    end

    return pade_table

end



pade_tab=newton_divided_difffunction(1,[1],2,[pi,3,-1],0,[1,2,5,1])
p,q=pade_tab[3]
test=polyval(p,2)/polyval(q,2)
##
q=newton_polynomial(1,[1],2,[0,3,-1],0,[1,2,5,1])





#TODO: write testing scheme
