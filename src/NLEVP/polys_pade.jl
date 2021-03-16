module Pade
import SpecialFunctions
## Polynomial type
"""
Polynomial type
"""
struct Polynomial #TODO: make type parametric use eltype(f)
 coeffs
 #constructor
    function Polynomial(coeffs)
        N=length(coeffs)
        n=0
        for i=N:-1:1
         if coeffs[i]==0
             n+=1
         else
             break
         end
        end
        return new(coeffs[1:(N-n)])
    end
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
        if p.coeffs[1]!=0
            txt="$(p.coeffs[1])"
        else
            txt=""
        end
        for n=2:N
            coeff=string(p.coeffs[n])
            if coeff==0
                continue
            end
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
    return +(a,-1*b)
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
#TODO: Devlop bigfloat prodrange
# function derive(a::Polynomial,k::Int=1)
#
#     N=length(a.coeffs)
#     if k>=N
#         return Polynomial([]) #TODO: type
#     end
#     d=zeros(typeof(a.coeffs[1]),N-k)
#     for i=1:N-k
#         d[i]=a.coeffs[i+k]*prodrange(i-1,i+k-1)
#     end
#     return Polynomial(d)
# end

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

# tests
# a=Polynomial([0.0,1.0])
# b=a*a
# c=derive(a)
# d=a+b+c
# e=d^3
# f=shift(e,2)
# e(2)-f(0)
# rational Polynomial
## Rational Polynomials
struct Rational_Polynomial
    num::Polynomial
    den::Polynomial

    #default constructor
    function Rational_Polynomial(num::Polynomial,den::Polynomial)
        num=copy(num.coeffs)
        den=copy(den.coeffs)
        if isempty(num) || num==[0]
            den=[1]
        elseif den[1]!=0
            den./=den[1]
            num./=den[1]
        end
        #TODO: else case
        return new(Polynomial(num),Polynomial(den))
    end
end
function (r::Rational_Polynomial)(z,k::Int=0)
    for i=1:k
        r=derive(r)
    end
    return r.num(z)/r.den(z)
end

function derive(a::Rational_Polynomial)
    return Rational_Polynomial(derive(a.num)*a.den-a.num*derive(a.den),a.den^2)
end
##



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
function derive(p::Polynomial)
    N=length(p.coeffs)
    return Polynomial([i*p.coeffs[i+1] for i=1:N-1] )
end

function pade(p::Polynomial,L,M) #TODO: harmonize symbols for eigenvalues, modes etc
  p=p.coeffs
  A=zeros(ComplexF64,M,M)
  for i=1:M
    for j=1:M
      if L+i-j>=0
        A[i,j]=p[L+i-j+1] #+1 is for 1-based indexing
      end
    end
  end
  b=A\(-p[L+2:L+M+1])
  b=[1;b]
  a=zeros(ComplexF64,L+1)
  for l=0:L
    for m=0:l
      if m<=M
        a[l+1]+=p[l-m+1]*b[m+1] #+1 is for zero-based indexing
      end
    end
  end
  a,b=Polynomial(a),Polynomial(b)
  return Rational_Polynomial(a,b)
end
end# module Pade
