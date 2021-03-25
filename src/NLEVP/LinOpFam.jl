## NLEVP
"""
    Term{T}

Single term of a linear operator family.

# Fields
- `coeff::T`: matrix coefficient of the term
- `func::Tuple`: tuple of scalar-valued functions multiplying the matrix coefficient
- `symbol::String`: character string for displaying the function(s)
- `params::Tuple`: tuple of tuples containing function symbols for each function
- `operator:String`: String for displaying the matrix coefficient
- `varlist::Array{Symbol,1}`: Array containing all symbols of fucntion arguments
that are used at least once
"""
struct Term{T}
  coeff::T
  func::Tuple
  symbol::String
  params::Tuple
  operator::String
  varlist::Array{Symbol,1} #TODO: this is redundant information
end

#constructor
function Term(coeff,func::Tuple,params::Tuple,symbol::String,operator::String)
  varlist=Symbol[]
  for par in params
    for var in par
      push!(varlist,var)
    end
  end
  unique!(varlist)
  return Term(coeff,func,symbol,params,operator, varlist)
end

"""
    term=Term(coeff,func::Tuple,params::Tuple,operator::String)

Standard constructor for type `Term`.

# Arguments
- `coeff`: matrix coefficient of the term
- `func::Tuple`: tuple of scalar-valued functions multiplying the matrix coefficient
- `params::Tuple`: tuple of tuples containing function symbols for each function
- `operator:String`: String for displaying the matrix coefficient

# Notes

The rendering of the functions of `term` ist automotized. If the passed functions
implement a method `f(z::Symbol)` then this method is used, otherwise the
functions will be simply displayed with the string `f`. You can overwrite this
behavior by explicitly passing a string for the function display using the syntax

  term =Term(coeff,func::Tuple,params::Tuple,symbol::String,operator::String)

where the extra argument `symbol`is the string used to display the function.

See also: [LinearOperatorFamily](@ref)
"""
function Term(coeff,func::Tuple,params::Tuple,operator::String)
  symbol=""
  for (f,p) in zip(func,params)
    if applicable(f,p...)
      symbol*=f(p...)
    else
      symbol*="f("
      for idx = 1:length(p)-1
        symbol*="$(p[idx]),"
      end
      symbol*="$(p[end]))"
    end
  end
  return Term(coeff,func,params,symbol,operator)
end


#Solution object
"""
    Solution

Type for the solution of an iterative eigensolver.

# Fields
- `params`: Dictionary mapping parameter symbols to their values
- `v`: right (direct) eigenvector
- `v_adj`: left (adjoint) eigenvector
- `eigval`: the symbol of the parameter that is the eigenvalue
- `eigval_pert`: extra data for asymptotic series expansion of the eigenvalue
- `v_pert`: extra data for asymptotic series expansion of the eigenvector
- `auxval`: the symbol of the parameter that is the auxiliary eigenvalue

See also: [LinearOperatorFamily](@ref)
"""
mutable struct Solution #IDEA: make immutable and parametrized type
  params
  v
  v_adj
  eigval
  eigval_pert
  v_pert
  auxval
end
#constructor
"""
    sol=Solution(params,v,v_adj,eigval,auxval=Symbol())

Standard constructor for a solution object. See [Solution](@ref) for details.
"""
function Solution(params,v,v_adj,eigval,auxval=Symbol())
  return Solution(deepcopy(params),v,v_adj,eigval,Dict{Symbol,Any}(),Dict{Symbol,Any}(),auxval) #TODO: copy params?
end

"""
    LinearOperatorFamily

Type for a linear operator family. The type is mutable.

# Fields
- `terms`: Array of terms forming the operator family
- `params`: Dictionary mapping parameter symbols to their values
- `eigval`: the symbol of the parameter that is the eigenvalue
- `auxval`: the symbol of the parameter that is the auxiliary eigenvalue
- `active`: list of symbols that can be actively change when using an instance of
the family as a function.
- `mode`: mode that is controlling which terms are evalauted when calling an
instance of the family. This field is not to be modified by the user.

See also: [Solution](@ref), [Term](@ref)
"""
mutable struct LinearOperatorFamily #TODO: add example to the doc
  terms::Array{Term,1}
  params::Dict{Symbol,ComplexF64}
  eigval::Symbol
  auxval::Symbol
  active::Array{Symbol,1}
  mode::Symbol
end

#constructors
"""
    L=LinearOperatorFamily(params,values)

Standard constructor for an empty Linear operator family.

# Arguments
- `params`: List of parameter symbols
- `values`: List of corresponding parameter values

# Note
The first parameter appearing in `params` will be assigned as eigenvalue. If
there are more than 1 parameters in the list, the last parameter will be designated
as auxiliary eigenvalue. (These choices can be changed after construction)
If `values`is omitted, then all parameters will be initialized with `NaN+NaN*im`.
If also `params` is omitted the family will be initialized with `:λ` as its
eigenvalue an no auxiliary eigenvalue.

Terms can be added to the family using the `+` operator or the more memory efficient
`push!` function. For instance `L+=term` or `push!(L,term)` both add `term`
to the list of terms.

See also: [Solution](@ref), [Term](@ref)
"""
# standard constructors
function LinearOperatorFamily(params,values)
  terms=Term[]
  eigval=Symbol(params[1])
  if length(params)>1
    auxval=Symbol(params[end])
  else
    auxval=Symbol("")
  end
  active=[eigval]
  pars=Dict{Symbol,ComplexF64}()
  for (p,v) in zip(params,values) #Implement alphabetical sorting
    pars[Symbol(p)]=v
  end
  mode=:all
  return LinearOperatorFamily(terms,pars,eigval,auxval,active,mode)
end
function LinearOperatorFamily()
  return LinearOperatorFamily(["λ",],[NaN+NaN*1im,])
end
function LinearOperatorFamily(params)
  return LinearOperatorFamily(params,[NaN+NaN*1im for a in params])
end

#loader
"""
    L=LinearOperatorFamily(fname::String)

Load and construct `LinearOperatorFamily` from file `fname`.

See also: [`save`](@ref)
"""
function LinearOperatorFamily(fname::String)
  D=read_toml(fname)
  eigval=D["eigval"]
  auxval=D["auxval"]
  params=[]
  vals=[]
  for (par,val) in D["params"]
    push!(params,par)
    push!(vals,val)
  end
  L=LinearOperatorFamily(params,vals)
  L.eigval=eigval
  L.auxval=auxval
  L.active=[eigval]
  for idx =1:length(D["/terms"])
    term=D["/terms"]["/"*string(idx)]
    I=term["/sparse_matrix"]["I"]
    J=term["/sparse_matrix"]["J"]
    V=term["/sparse_matrix"]["V"]
    m, n = term["size"]
    M=sparse(I,J,V,m,n)
    push!(L,Term(M,term["functions"],term["params"],term["symbol"],term["operator"]))
  end
  return L
end

import Dates
#saver
"""
    save(fname::String,L::LinearOperatorFamily)

Save `L` to file `fname`. The file is utf8 encoded and adheres to a Julia-enriched TOML standard.

See also: [`LinearOperatorFamily`](@ref)
"""
function save(fname,L::LinearOperatorFamily)
  date=string(Dates.now(Dates.UTC))

  eq=""
  for term in L.terms
    if term.operator[1]=='_' #TODO: put this into signature?
      continue
    end
    eq*="+"*string(term)
  end

  open(fname,"w") do f
    write(f,"# LinearOperatorFamily version 0\n")
    write(f,"#"*date*"\n")
    write(f,"#"*eq*"\n")
    write(f,"params=[")
    for (key,value) in L.params
        write(f,"(:$(string(key)),$value),\n")
    end
    write(f,"]\n")

    write(f,"eigval=:$(L.eigval)\n")
    write(f,"auxval=:$(L.auxval)\n")
    #TODO activitiy and terms
    write(f,"[terms]\n")
    for (idx,term) in enumerate(L.terms)
      write(f,"\t[terms.$idx]\n")
      write(f,"\tfunctions=(")
      for func in term.func
        write(f,"$func,")
      end
      write(f,")\n")
      write(f,"\tsymbol=\"$(term.symbol)\"\n")
      write(f,"\tparams=$(term.params)\n")
      # for param in term.params
      #   write(f,"[")
      #   for par in param
      #     write(f,":$(String(par)),")
      #   end
      #   write(f,"],")
      # end
      # write(f,"]\n")
      write(f,"\toperator=\"$(term.operator)\"\n")
      m,n=size(term.coeff)
      write(f,"\tsize=[$m,$n]\n")
      #write(f,"\tis_sparse=$(SparseArrays.issparse(term.coeff))\n")
      write(f,"\t\t[terms.$idx.sparse_matrix]\n")
      I,J,V=SparseArrays.findnz(term.coeff)
      write(f,"\t\tI=$I\n")
      write(f,"\t\tJ=$J\n")
      #write(f,"\t\tV=$V\n")
      write(f,"\t\tV=Complex{Float64}[")
      for v in V
        if imag(v)>=0
          write(f,"$(real(v))+$(imag(v))im,")
        else
          write(f,"$(real(v))$(imag(v))im,")
        end
      end
      write(f,"]\n")
      write(f,"\n")
    end
  end
end

#TODO lesen
#eval funktion nutzen
#beginnt ist mit :  ????
#beginnt es mit (
#beginnt es mit [   ???
#dann eval
# es gibt nur drei arten: tags, variablen, listen/tuple rest ist eval

import Base.push!
function push!(L::LinearOperatorFamily,T::Term)
  d=Dict()
  for (idx,term) in enumerate(L.terms)
    d[(term.func, term.params)]=idx
  end
  signature=(T.func,T.params)
  if signature in keys(d) #change existing term if signature known
    idx=d[signature]
    coeff=L.terms[idx].coeff+T.coeff
    if LinearAlgebra.norm(coeff)==0
      deleteat!(L.terms,idx) #delte if resulting coeff is 0
      #check for unbound variables and delete
      vars=[]
      for term in L.terms
        for pars in term.params
          for par in pars
            if par ∉ vars
              push!(vars,par)
            end
          end
        end
      end
      for par in keys(T.params)
        if par ∉ vars
          delete!(L.params,par)
        end
      end

    else
      L.terms[idx]=Term(coeff,L.terms[idx].func,L.terms[idx].symbol,L.terms[idx].params,L.terms[idx].operator,L.terms[idx].varlist) #overwrite term if coeff is non-zero
    end
  else #add term if signature is new
    for pars in T.params
      for par in pars
        if par ∉ keys(L.params)
          L.params[par]=NaN+NaN*1im #initialize variable if its new
        end
      end
    end
    push!(L.terms,T)
  end
end
# function push!(a::LinearOperatorFamily,b::LinearOperatorFamily)
#    push!(a.terms,b.terms)
# end
#
#
import Base.(+)
function (+)(a::LinearOperatorFamily,b::Term)
  L=deepcopy(a)
  push!(L,b)
  return L
end
function (+)(b::Term,a::LinearOperatorFamily,)
  L=deepcopy(a)
  push!(L,b)
  return L
end
import Base.(-)
function (-)(a::LinearOperatorFamily,b::Term)
  L=deepcopy(a)
  push!(L,Term(-b.coeff,b.func,b.symbol,b.params,b.operator,b.varlist))
  return L
end
function (-)(b::Term,a::LinearOperatorFamily)
  L=deepcopy(a)
  for i=1:length(L.terms)
    L.terms[i].coeff*=-1
  end
  push!(L,b)
  return L
end




#nice display of LinearOperatorFamilies in interactive and other modes
import Base.show
function show(io::IO,L::LinearOperatorFamily)
  if !isempty(L.terms)
    shape=size(L.terms[1].coeff)
    if length(shape)==2
      txt="$(shape[1])×$(shape[2])-dimensional operator family: \n\n"
    else
      txt="$(shape[1])-dimensional vector family: \n\n"
    end
  else
    txt="empty operator family\n\n"
  end

  eq=""
  for term in L.terms
    if term.operator[1]=='_'
      continue
    end
    eq*="+"*string(term)
  end

  parameter_list="\n\nParameters\n----------\n"
  for (par,val) in L.params
    parameter_list*=string(par)*"\t"*string(val)*"\n"
  end
  print(io, txt*eq[2:end]*parameter_list)
end

import Base.string
function string(T::Term)
  if T.symbol==""
    txt=""
  else
    txt=string(T.symbol)*"*"
  end
  txt*=T.operator
end
function show(io::IO,T::Term)
  print(io,string(T))
end

#make solution showable
function string(sol::Solution)
  txt="""####Solution####
  eigval:
  $(sol.eigval) = $(sol.params[sol.eigval])

  Parameters:
  """
  for (key,val) in sol.params
    if key!=sol.eigval && key!=sol.auxval
      txt*="$key = $val\n"
    end
  end

  if sol.auxval in keys(sol.params)
    txt*="""

    Residual:
    abs($(sol.auxval)) = $(abs(sol.params[sol.auxval]))
    """
  end
  return txt
end
function show(io::IO,sol::Solution)
  print(io,string(sol))
end


#make terms callable
function (term::Term)(dict::Dict{Symbol,Tuple{ComplexF64,Int64}})
  coeff::ComplexF64=1 #TODO parametrize type
  for (func,pars) in zip(term.func, term.params)
    args=[]
    dargs=[]
    for par in pars
      push!(args,dict[par][1])
      push!(dargs,dict[par][2])
    end

    coeff*=func(args...,dargs...)
  end
  return coeff*term.coeff
end

#make LinearOperatorFamily callable
function(L::LinearOperatorFamily)(args...;oplist=[],in_or_ex=false)
  if L.mode==:all # if mode is all the first n args correspond to the values of the n active variables
    for (var,val) in zip(L.active,args)
      L.params[var]=val #change the active variables
    end
  end
  if L.mode==:all && length(args)==length(L.active)
    derivs=zeros(Int64,length(L.active))
  else
    derivs=args[end-length(L.active)+1:end] #TODO implement sanity check on length of args
  end

  derivDict=Dict{Symbol,Int64}()
  for (var, drv) in zip(L.active,derivs)
    derivDict[var]=drv #create a dictionairy for the active variables
  end

  coeff=spzeros(size(L.terms[1].coeff)...)  #TODO: improve this to handle more flexible matrices

  for term in L.terms
    if (!in_or_ex && term.operator in oplist) || (in_or_ex && !(term.operator in oplist)) || (L.mode!=:householder && term.operator=="__aux__")  #TODO: consider deprecating this feature together with nicoud and picard
      continue
    end

    #check whether term is constant w.r.t. to some parameter then deriv is 0 thus continue
    skip=false
    for (var,d) in zip(L.active, derivs)
      if d>0 && !(var in term.varlist)
        skip=true
        break
      end
    end
    if skip
      continue
    end
    dict=Dict{Symbol,Tuple{Complex{Float64},Int64}}()
    for var in term.varlist
      dict[var]=(L.params[var], var in keys(derivDict) ? derivDict[var] : 0)
    end
    coeff+=term(dict) #
  end

  if L.mode in [:compact,:householder]
    coeff/=prod(factorial.(float.(args[end-length(L.active)+1:end])))
  end

  return coeff
end

#wrapper to perturbation theory
#TODO: functional programming renders this obsolete, no?
"""
    perturb!(sol::Solution,L::LinearOperatorFamily,param::Symbol,N::Int; <keyword arguments>)

Compute the `N`th order power series perturbation coefficients for the solution `sol` of the nonlineaer eigenvalue problem given by the operator family `L` with respect to the parameter `param`. The coefficients will be stored in the field `sol.eigval_pert` and `sol.v_pert` for the eigenvalue and the eignevector, respectively.

# Keyword Arguments
- `mode = :compact`: parameter controlling internal programm flow. Use the default, unless you know what you are doing.

# Notes
For large perturbation orders `N` the method might be slow.

See also: [`perturb_fast!`](@ref)
"""
 function perturb!(sol::Solution,L::LinearOperatorFamily,param::Symbol,N::Int;mode=:compact)
   active=L.active #TODO: copy?
   params=L.params
   L.params=sol.params
   current_mode=L.mode
   L.active=[sol.eigval, param]
   L.mode=mode
   pade_symbol=Symbol("$(string(param))/Taylor")
   sol.eigval_pert[pade_symbol],sol.v_pert[pade_symbol]=perturb(L,N,sol.v,sol.v_adj)
   sol.eigval_pert[pade_symbol][1]=sol.params[sol.eigval]
   L.active=active
   L.mode=current_mode
   L.params=params
   return
 end

"""
    perturb_fast!(sol::Solution,L::LinearOperatorFamily,param::Symbol,N::Int; <keyword arguments>)

Compute the `N`th order power series perturbation coefficients for the solution `sol` of the nonlineaer eigenvalue problem given by the operator family `L` with respect to the parameter `param`. The coefficients will be stored in the field `sol.eigval_pert` and `sol.v_pert` for the eigenvalue and the eignevector, respectively.

# Keyword Arguments
- `mode = :compact`: parameter controlling internal programm flow. Use the default, unless you know what you are doing.

# Notes
This method reads multi-indeces for the computation of the power series coefficients from disk. Make sure that JulHoltz is properly installed to use this fast method.

See also: [`perturb!`](@ref)
"""
 function perturb_fast!(sol::Solution,L::LinearOperatorFamily,param::Symbol,N::Int;mode=:compact)
   active=L.active #TODO: copy?
   params=L.params
   L.params=sol.params
   current_mode=L.mode
   L.active=[sol.eigval, param]
   L.mode=mode
   pade_symbol=Symbol("$(string(param))/Taylor")
   sol.eigval_pert[pade_symbol],sol.v_pert[pade_symbol]=perturb_disk(L,N,sol.v,sol.v_adj)
   sol.eigval_pert[pade_symbol][1]=sol.params[sol.eigval]
   L.active=active
   L.mode=current_mode
   L.params=params
   return
 end

"""
    perturb_norm!(sol::Solution,L::LinearOperatorFamily,param::Symbol,N::Int; <keyword arguments>)

 Compute the `N`th order power series perturbation coefficients for the solution `sol` of the nonlineaer eigenvalue problem given by the operator family `L` with respect to the parameter `param`. The coefficients will be stored in the field `sol.eigval_pert` and `sol.v_pert` for the eigenvalue and the eignevector, respectively.

 # Keyword Arguments
 - `mode = :compact`: parameter controlling internal programm flow. Use the default, unless you know what you are doing.

 # Notes
 This method reads multi-indeces for the computation of the power series coefficients from disk. Make sure that JulHoltz is properly installed to use this fast method.

 See also: [`perturb!`](@ref)
 """
  function perturb_norm!(sol::Solution,L::LinearOperatorFamily,param::Symbol,N::Int;mode=:compact)
    active=L.active #TODO: copy?
    params=L.params
    L.params=sol.params
    current_mode=L.mode
    L.active=[sol.eigval, param]
    L.mode=mode
    pade_symbol=Symbol("$(string(param))/Taylor")
    sol.eigval_pert[pade_symbol],sol.v_pert[pade_symbol]=perturb_norm(L,N,sol.v,sol.v_adj)
    sol.eigval_pert[pade_symbol][1]=sol.params[sol.eigval]
    L.active=active
    L.mode=current_mode
    L.params=params
    return
  end

 #TODO: Implement solution class

function pade(ω,L,M) #TODO: harmonize symbols for eigenvalues, modes etc
  A=zeros(ComplexF64,M,M)
  for i=1:M
    for j=1:M
      if L+i-j>=0
        A[i,j]=ω[L+i-j+1] #+1 is for 1-based indexing
      end
    end
  end
  b=A\(-ω[L+2:L+M+1])
  b=[1;b]
  a=zeros(ComplexF64,L+1)
  for l=0:L
    for m=0:l
      if m<=M
        a[l+1]+=ω[l-m+1]*b[m+1] #+1 is for zero-based indexing
      end
    end
  end
  return a,b
end



function pade!(sol::Solution,param,L::Int64,M::Int64;vector=false)
  #TODO: implement sanity check whether L+M+1=length(sol.eigval_pert[:Taylor])
  pade_symbol=Symbol("$(string(param))/[$L/$M]")
  taylor_symbol=Symbol("$(string(param))/Taylor")
  sol.eigval_pert[pade_symbol]=pade(sol.eigval_pert[taylor_symbol],L,M)

  #TODO Degeneracy
  if vector
    D=length(sol.v)
    #preallocation
    A=Array{Array{ComplexF64}}(undef,L+1)
    B=Array{Array{ComplexF64}}(undef,M+1)
    for idx in 1:length(A)
      A[idx]=Array{ComplexF64}(undef,D)
    end
    for idx in 1:length(B)
      B[idx]=Array{ComplexF64}(undef,D)
    end
    for idx in 1:D
      v=Array{ComplexF64}(undef,L+M+1)
      for ord = 1:L+M+1
        v[ord]=sol.v_pert[taylor_symbol][ord][idx]
      end
      a,b=pade(v,L,M)
      for ord=1:length(A)
        A[ord][idx]=a[ord]
      end
      for ord=1:length(B)
        B[ord][idx]=b[ord]
      end
    end
    sol.v_pert[pade_symbol]=A,B
  end
  return
end


#make solution object callable
function (sol::Solution)(param::Symbol,ε,L=0,M=0; vector=false) #Todo not really performant, consider syntax that allows for vectorization in eigenvector
    pade_symbol=Symbol("$(string(param))/[$L/$M]")
    if pade_symbol ∉ keys(sol.eigval_pert) || (vector && pade_symbol ∉ keys(sol.v_pert))
      pade!(sol,param,L,M,vector=vector)
    end
    a,b=sol.eigval_pert[pade_symbol]
    Δε=ε-sol.params[param]
    eigval=polyval(a,Δε)/polyval(b,Δε)
    if !vector
      return eigval
    else
      A, B = sol.v_pert[pade_symbol]
      eigvec = polyval(A,Δε) ./ polyval(B,Δε)
      return eigval, eigvec
    end
end
#TODO: implement parameter activity in solution type

# function polyval(A,z)
#     f=zeros(eltype(A[1]),length(A[1]))
#     for (idx,a) in enumerate(A)
#         f.+=a*(z^(idx-1))
#     end
#     return f
# end

#polyval(p,z)= sum(p.*(z.^(0:length(p)-1)))  #TODO: Although, this is not performance critical. Consider some performance aspects
"""
    f=polyval(p,z)
Evaluate polynomial f(z)=∑_i p[i]z^1 at z, using Horner's scheme.
"""
function polyval(p,z)
  f=ones(size(z))*p[end]
  for i = length(p)-1:-1:1
     f.*=z
     f.+=p[i]
  end
  return f
end
function polyval(p,z::Number)
  f=p[end]
  for i = length(p)-1:-1:1
     f*=z
     f+=p[i]
  end
  return f
end





function estimate_pol(ω::Array{Complex{Float64},1})
  N=length(ω)
  Δε=zeros(ComplexF64,N-2)
  k=zeros(ComplexF64,N-2)
  for j =2:length(ω)-1
    i=j-1 # index shift to account for 1-based indexing
    denom=(i+1)*ω[j+1]*ω[j-1]-i*ω[j]^2
    Δε[i]=ω[j]*ω[j-1]/denom
    k[i]=(i^2-1)*ω[j+1]*ω[j-1]-(i*ω[j])^2
  end
  return Δε,k
end

function estimate_pol(sol::Solution,param::Symbol)
  pade_symbol=Symbol("$(string(param))/Taylor")
  return estimate_pol(sol.eigval_pert[pade_symbol])
end

function conv_radius(a::Array{Complex{Float64},1})
  N=length(a)
  r=zeros(Float64,N-1)
  for n=1:N-1
    r[n]=abs(a[n]/a[n+1])
  end
  return r
end

function conv_radius(sol::Solution,param::Symbol)
  pade_symbol=Symbol("$(string(param))/Taylor")
  return conv_radius(sol.eigval_pert[pade_symbol])
end
