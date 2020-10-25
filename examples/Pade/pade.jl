#This file will be turned in more flexible algebra routines

##header
import SpecialFunctions

## polynomial type
"""
polynomial type
"""
 struct polynomial #TODO: make type parametric use eltype(f)
     coeffs
 end
#make polynomial callable with Horner-scheme evaluation
function (p::polynomial)(z)
    f=z*p.coeffs[end]
    for i = length(p.coeffs)-1:-1:1
       f*=z
       f+=p.coeffs[i]
    end
    return f
end

import Base.+
import Base.-
import Base.*
import Base.^



function +(a::polynomial,b::polynomial)
    if length(a.coeffs)>length(b.coeffs)
        a,b=b,a
    end
    c=copy(b.coeffs)
    for (idx,val) in enumerate(a.coeffs)
        c[idx]+=val
    end
    return polynomial(c)
end

function -(a::polynomial,b::polynomial)
    return +(a,-b)
end

function *(a::polynomial,b::polynomial)
    I=length(a.coeffs)
    J=length(b.coeffs)
    K=I+J-1
    c=zeros(typeof(a.coeffs[1]),K) #TODO: sanity check for type of b
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
    return polynomial(c)
end

function ^(p::polynomial, k::Int)
    b=polynomial(copy(p.coeffs))
    for i=1:k-1
        b*=p
    end
    return b
end

function prodrange(start,stop)
    return SpecialFunctions.factorial(stop)/SpecialFunctions.factorial(start)
end

"""
    b=derive(a::polynomial,k::Int=1)

Compute `k`th derivative of polynomial `a`.
"""
function derive(a::polynomial,k::Int=1)
    N=length(a.coeffs)
    d=zeros(typeof(a.coeffs[1]),N-k)
    for i=1:N-k
        d[i]=a.coeffs[i+k]*prodrange(i-1,i+k-1)
    end
    return polynomial(d)
end

"""
    g=shift(f,Δz)

Shift polynomial with respect to its argument. (Taylor shift)
"""
function shift(f,Δz)
    g=Array{eltype(f)}(undef,length(f.coeffs))
    for n=0:length(f.coeffs)-1
        g[n+1]=f(Δz)/SpecialFunctions.factorial(n*1.0)
        f=polynomial(f.coeffs[2:end])
        f.coeffs.*=1:length(f.coeffs)
    end
    return polynomial(g)
end

## tests
a=polynomial([0.0,1.0])
b=a*a
c=derive(a)
d=a+b+c
e=d^3
f=shift(e,2)
e(2)-f(0)
## rational polynomial
struct rational_polynomial
    num::polynomial
    den::polynomial
end

function derive(a::rational_polynomial)
    return rational_polynomial(derive(a.num)*a.den-a.num*derive(a.den),a.den^2)
end

"""
    derive(a::Any, k:Int)

Take the `k`th derivative of `a` by calling the derive function multiple times.

# Notes
This is a fallback implementation.
"""
function derive(a::Any,k:Int)
    for i=1:k
        a=derive(a)
    end
    return a
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
    basis_function=polynomial([1])
    idx=comb[1:1]


    c=polynomial(coeffs[idx])
    p=basis_function*c
    for n=2:length(comb)
        idx=comb[1:n]
        c=polynomial(coeffs[idx])
        basis_function=basis_function*polynomial(-z[idx[end-1]],1]
        p=p+ᴾc*ᴾbasis_function
    end
    basis_function=basis_function*ᴾ[-z[comb[end]],1]


    #Kronecker's Algorithm to find the interpolating Pade approximant
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
