## algebra
function pow0(z::ComplexF64,k::Int=0)::ComplexF64
  if k==0
    return 1
  elseif k>0
    return 0
  else
    return complex(NaN)
  end
end

function pow1(z::ComplexF64,k::Int=0)::ComplexF64
  if k==0
    return z
  elseif k==1
    return 1
  elseif k>1
    return 0
  else
    return complex(NaN)
  end
end

function pow2(z::ComplexF64,k::Int=0)::ComplexF64
  if k==0
    return z^2
  elseif k==1
    return 2*z
  elseif k==2
    return 2
  elseif k>2
    return 0
  else
    return complex(NaN)
  end
end

#TODO turn this into a macro
function pow(z::ComplexF64,k::Int64,a::Int64)::ComplexF64
  if k>a
    f= 0
  elseif k>=0
    f=1
    i=a
    for j=0:k-1
      f*=i
      i-=1
    end
    f*=z^(a-k)
  else
    f=complex(NaN)
  end
  return f
end

function pow_a(a::Int64)
  return (z::ComplexF64,k::Int64=0)-> pow(z,k,a)
end

function exp_delay(ω,τ,m,n)
  a=-1.0im
  u(z,l)=pow(z,l,m)
  f=0
  for i = 0:n
    f+=binomial(n,i)*u(τ,i)*(a*ω)^(n-i)
  end
  f*=a^m*exp(a*ω*τ)
  return f
end

tau_delay=exp_delay

#TODO create a macro for generating functions

function exp_az(z::ComplexF64,a::ComplexF64,k::Int)::ComplexF64
  if k>=0
    return a^k*exp(a*z)
  else
    return nothing
  end
end


# function powa(z::ComplexF64,k::Int,a::Int)
#   if
# end
function z_exp_iaz(z::ComplexF64, a::ComplexF64,m=0,n=0)::ComplexF64
  if m==n==0
    return z*exp(1im*a*z)
  elseif m==1
    (1im*a*z +1)*exp(1im*a*z)
  elseif n==1
    return 1im*z^2*exp(1im*a*z)
  else
    print("Argh!!!!!")
  end
end

function z_exp__iaz(z::ComplexF64, a::ComplexF64,m=0,n=0)::ComplexF64
  if m==n==0
    return z*exp(-1im*a*z)
  elseif m==1
    (-1im*a*z +1)*exp(-1im*a*z)
  elseif n==1
    return- 1im*z^2*exp(-1im*a*z)
  else
    print("Argh!!!!!")
  end
end

function exp_pm(s)
a=s*1.0im
  function exp_delay(ω,τ,m,n)
    u(z,l)=pow(z,l,m)
    f=0
    for i = 0:n
      f+=binomial(n,i)*u(τ,i)*(a*ω)^(n-i)
    end
    f*=a^m*exp(a*ω*τ)
    return f
  end
return exp_delay
end

function exp_ax2(z,a,n)
  if a==0.0+0im
    if n==0
      return 1.0+0.0im
    else
      return 0.0+0.0im
    end
  end

  f=0.0+0.0im
  A=a^n
  Z=z^n
  cnst=2^n*factorial(n)
  for k=0:n÷2
    coeff=0.0+0.0im
    coeff+=cnst*4.0^(-k)/factorial(k)/factorial(n-2*k)
    coeff*=A
    coeff*=Z
    f+=coeff
    A/=a
    Z/=z^2
  end
  f*=exp(a*z^2)
  return f
end

function exp_ax2mxit(z,τ,a,m,n,k)
  f=pow_a(n+2*k)
  g(z,l)=exp_ax2(z,a,l)
  h(z,l)=exp_delay(z,τ,l,0)
  coeff=0.0+0.0im
  multinomialcoeff=factorial(m)
  for ii=0:m
    multi_ii=factorial(ii)
    coeff_ii=h(z,ii)
    for jj=0:m-ii
      multi_jj=factorial(jj)
      coeff_jj=g(z,jj)
      kk=m-jj-ii
      multi=multinomialcoeff/multi_ii/multi_jj/factorial(kk)
      coeff+=multi*f(z,kk)*coeff_jj*coeff_ii
    end
  end
    coeff*=(-1.0im)^n
    return coeff
end
