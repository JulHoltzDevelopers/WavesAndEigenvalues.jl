##Ok let's go and write an iterator for Kelleher's algorithm
struct part
  n
end

import Base
#import IterativeSolvers
import SpecialFunctions

#initializer
function Base.iterate(iter::part)
  a=zeros(Int64,iter.n)
  k=2
  y=iter.n-1
  if k!=1
    x=a[k-1]+1
    k-=1
    while 2*x <=y
      a[k]=x
      y-=x
      k+=1
    end
    l=k+1
    if x<=y
      a[k]=x
      a[l]=y
      #put!(c,a[1:k+1])
      x+=1
      y-=1
      return a[1:k+1],(a,k,l,x,y)
    end
    a[k]=x+y
    y=x+y-1
    ## header again there is probably an easier way to code this
    p= a[1:k]
    if k!=1
      x=a[k-1]+1
      k-=1
      while 2*x <=y
        a[k]=x
        y-=x
        k+=1
      end
      l=k+1
    else
      k=-1 # sentinel value see below
    end
    return p,(a,k,l,x,y)
  end
end

function Base.iterate(iter::part,state)
  a,k,l,x,y=state
  if k==-1; return nothing; end #minus 1 is a sentinel value
  if x<=y
    a[k]=x
    a[l]=y
    #put!(c,a[1:k+1])
    x+=1
    y-=1
    return a[1:k+1],(a,k,l,x,y)
  else
    a[k]=x+y
    y=x+y-1
    p=a[1:k]# copy next output
    if k!=1
      x=a[k-1]+1
      k-=1
      while 2*x <=y
        a[k]=x
        y-=x
        k+=1
      end
      l=k+1
    else
      k=-1 #set sentinel value for termination
    end
    return p,(a,k,l,x,y)
  end
end

# #small performance test
# #Watch this PyHoltz
# zeit=time()
# for (num,p) in enumerate(part(15))
#   println(num,p)
# end
# zeit=time()-zeit
# print("iterator:",zeit)





function part2mult(p)
  z=sum(p)
  mu=zeros(Int64,z)
  if p !=[0]
    for i in p
      mu[i]+=1
    end
  end
  return mu
end

# #next tiny test
# for (num,p) in enumerate(part(5))
#   println(num,p,part2mult(p))
# end

function multinomcoeff(mu)
  return SpecialFunctions.factorial(float(sum(mu)))/prod(SpecialFunctions.factorial.(float.(mu)))
end

function weigh(mu)
  weight=0
  for (g,mu_g) in enumerate(mu)
    weight+=g*mu_g #TODO: check type
  end
  return weight
end



function generate_MN(N)
  #MN=[Dict{Tuple{Int64,Int64},Array{Array{Int64,1},1}}()]
  MN=[Dict{Array{Int64,1},Array{Array{Int64,1},1}}()]
  #MN=[Dict()]
  for k=1:N
    #MU=Dict{Tuple{Int64,Int64},Array{Array{Int64,1},1}}()
    MU=Dict{Array{Int64,1},Array{Array{Int64,1},1}}()
    #MU=Dict()
    for n=1:k
      key=[0,n]
      mu=[0]
      mu=part2mult(mu)
      #mu=Array{Int64,1}[]
      if key in keys(MU)
        push!(MU[key],mu)
      else
        MU[key]=[mu]
      end
    end

    for m = 1:k
      for mu in part(m)
        if mu == [k]
          continue
        end
        mu=part2mult(mu)
        for n=0:k-m
          key=[sum(mu),n]
          if key in keys(MU)
            push!(MU[key],mu)
          else
            MU[key]=[mu]
          end
        end
      end
    end
    push!(MN,MU)
  end
  return MN
end




function tuple2idx(tpl,ord)
  m,n=tpl
  return sum(ord-m+1:ord)+n+m
end

function mult2str(mu)
  txt=""
  for mu_g in mu
    txt*=string(mu_g)*" "
  end
  if length(txt)!=0
    txt=txt[1:end-1]
    txt*="\n"
  end
 return txt
end

function generate_multi_indeces_at_order(k;to_disk=false,compressed=false)
  if to_disk
    pack="../src/compressed_perturbation_data/" #This file location is relative to the deps directory as this is the one where build pkg is run.
    dir="$k/"
  else
    Mu=[ Array{Int16,1}[] for i in 1:tuple2idx((k,0),k) ] #TODO specify type
  end

  for n=1:k

    if to_disk && !ispath(pack*dir)
      mkpath(pack*dir)
    end
    key=(0,n)
    p=[0]
    mu=part2mult(p)
    out=compressed ? p :  mu
    if to_disk
      fname="$(key[1])_$(key[2])"
      open(pack*dir*fname,"a") do file
        write(file,mult2str(out))
      end
    else
      idx=tuple2idx(key,k)
      push!(Mu[idx],out)
    end
  end

  for m = 1:k
    for p in part(m)
      if p == [k]
        continue
      end
      mu=part2mult(p)
      out=compressed ? p :  mu
      for n=0:k-m
        key=(sum(mu),n)

        if to_disk
          fname="$(key[1])_$(key[2])"
          open(pack*dir*fname,"a") do file
            write(file,mult2str(out))
          end
        else
          idx=tuple2idx(key,k)
          push!(Mu[idx],out)
        end

      end
    end
  end

  if !to_disk
    return Mu
  end
end



function efficient_MN(N)
  MN=Array{Array{Array{Array{Int16,1},1},1}}(undef,N)
  #MN=[]
  for k =1:N
  #push!(MN,generate_multi_indeces_at_order(k))
  MN[k]=generate_multi_indeces_at_order(k)
  end
  return MN
end

function disk_MN(N)
  for k=1:N
    generate_multi_indeces_at_order(k,to_disk=true,compressed=true)
  end
end
##
# zeit=time()
# MN=generate_MN(50);
# zeit=time()-zeit
# println("zeit: ", zeit)
##
# function perturb(L,N,v0,v0Adj)
#   #normalize
#   v0/=sqrt(v0'*v0)
#   v0Adj/=v0Adj'*L(1,0)*v0
#   λ=Array{ComplexF64}(undef,N+1)
#   v=Array{Array{Complex{Float64}},1}(undef,N+1)
#   v[1]=v0
#   dim=size(v0)[1]
#   L00=L(0,0)
#   sparse = SparseArrays.issparse(L00)
#   L00=SparseArrays.lu(L(0,0),check=false) #L(0,0) is singular but in most of the cases we are lucky and a LU factorization exists
#   if sparse && L00.status!=0 #lu failed... TODO: implement LU check for dense arrays
#      L00=SparseArrays.qr(L(0,0)) #we try a QR factorization
#   end
#   #TODO consider using qr for all cases
#
#
#   for k=1:N
#     r=zeros(ComplexF64,dim)
#     for n = 1:k
#       r+=L(0,n)*v[k-n+1] #+1 to start indexing with zero
#     end
#     for m =1:k
#       for mu in part(m)
#         if mu ==[k]
#           continue
#         end
#         mu=part2mult(mu)
#         for n=0:k-m
#           coeff=1
#           for (g, mu_g) in enumerate(mu)
#             coeff*=λ[g+1]^mu_g #+1 because indexing starts with zero
#           end
#           r+=L(sum(mu),n)*v[k-n-m+1]*multinomcoeff(mu)*coeff
#         end
#       end
#     end
#     λ[k+1]=-v0Adj'*r/(v0Adj'*L(1,0)*v0)
#     #v[k+1]= IterativeSolvers.gmres(L(0,0),-(r+λ[k+1]*L(1,0)*v0)) This is crap!!!
#     #v[k+1]=L(0,0)\-(r+λ[k+1]*L(1,0)*v0) #This works but LU factorization should be done outside of the loop
#     v[k+1]=L00\-(r+λ[k+1]*L(1,0)*v0) #yeah !
#     v[k+1]-=(v0'*v[k+1])*v0 #if v[k+1] is orthogonal to v0 it features minimum norm
#     #
#     #println("status: $(L00.status)")
#     #res=L(0,0)*v[k+1]+(r+λ[k+1]*L(1,0)*v0)
#     #res=LinearAlgebra.norm(res)
#   end
#   return λ,v
# end

function perturb(L,N,v0,v0Adj)
  #normalize
  v0/=sqrt(v0'*v0)
  v0Adj/=v0Adj'*L(1,0)*v0
  λ=Array{ComplexF64}(undef,N+1)
  v=Array{Array{Complex{Float64}},1}(undef,N+1)
  v[1]=v0
  dim=size(v0)[1]
  L00=L(0,0)
  sparse = SparseArrays.issparse(L00)
  L00=SparseArrays.lu(L(0,0),check=false) #L(0,0) is singular but in most of the cases we are lucky and a LU factorization exists
  if sparse && L00.status!=0 #lu failed... TODO: implement LU check for dense arrays
     L00=SparseArrays.qr(L(0,0)) #we try a QR factorization
  end
  #TODO consider using qr for all cases


  for k=1:N
    r=zeros(ComplexF64,dim)
    for n = 1:k
      r+=L(0,n)*v[k-n+1] #+1 to start indexing with zero
    end
    for m =1:k
      for mu in part(m)
        if mu ==[k]
          continue
        end
        mu=part2mult(mu)
        for n=0:k-m
          coeff=1
          for (g, mu_g) in enumerate(mu)
            coeff*=λ[g+1]^mu_g #+1 because indexing starts with zero
          end
          r+=L(sum(mu),n)*v[k-n-m+1]*multinomcoeff(mu)*coeff
        end
      end
    end
    λ[k+1]=-v0Adj'*r/(v0Adj'*L(1,0)*v0)
    #v[k+1]= IterativeSolvers.gmres(L(0,0),-(r+λ[k+1]*L(1,0)*v0)) This is crap!!!
    #v[k+1]=L(0,0)\-(r+λ[k+1]*L(1,0)*v0) #This works but LU factorization should be done outside of the loop
    v[k+1]=L00\-(r+λ[k+1]*L(1,0)*v0) #yeah !
    v[k+1]-=(v0'*v[k+1])*v0 #if v[k+1] is orthogonal to v0 it features minimum norm
    #
    #println("status: $(L00.status)")
    #res=L(0,0)*v[k+1]+(r+λ[k+1]*L(1,0)*v0)
    #res=LinearAlgebra.norm(res)
  end
  return λ,v
end





#TODO: consider stronger compression
function perturb_disk(L,N,v0,v0Adj)
  #normalize
  zeit=time()
  v0/=sqrt(v0'*v0)
  v0Adj/=v0Adj'*L(1,0)*v0
  λ=Array{ComplexF64}(undef,N+1)
  v=Array{Array{Complex{Float64}},1}(undef,N+1)
  v[1]=v0
  dim=size(v0)[1]
  L00=L(0,0)
  sparse = SparseArrays.issparse(L00)
  L00=SparseArrays.lu(L(0,0),check=false) #L(0,0) is singular but in most of the cases were are lucky and a LU factorization exists
  if sparse && L00.status!=0 #lu failed... TODO: implement LU check for dense arrays
     L00=SparseArrays.qr(L(0,0)) #we try a QR factorization
  end
  #TODO consider using qr for all cases
  pack=(@__DIR__)*"/compressed_perturbation_data/"
  for k=1:N
    dir="$k/"
    r=zeros(ComplexF64,dim)
    for m=0:k
      for n=0:k-m
        if m==n==0 || k==m==1
          continue
        end
        fname="$(m)_$(n)"
        w=zeros(ComplexF64,dim)
        open(pack*dir*fname,"r") do file
          for str in eachline(file)
            #str=readline(file)
            mu=part2mult(parse.(Int64,split(str)))

            coeff=1
            for (g,mu_g) in enumerate(mu)
              coeff*=λ[g+1]^mu_g #+1 because indexing starts with zero
            end
            w+=v[k-n-weigh(mu)+1]*multinomcoeff(mu)*coeff
            end
        end
        r+=L(m,n)*w
      end
    end
    #println("##########")
    #println("r: $r")
    #v[k+1]= IterativeSolvers.gmres(L(0,0),-(r+λ[k+1]*L(1,0)*v0)) This is crap!!!
    λ[k+1]=-v0Adj'*r/(v0Adj'*L(1,0)*v0)
    #println("λ: $(λ[k+1])")
    #v[k+1]=L(0,0)\-(r+λ[k+1]*L(1,0)*v0) #This works but LU factorization should be done outside of the loop
    #println("rhs: $(r+λ[k+1]*L(1,0)*v0)")
    v[k+1]=L00\-(r+λ[k+1]*L(1,0)*v0) #yeah !
    v[k+1]-=(v0'*v[k+1])*v0 #if v[k+1] is orthogonal to v0 it features minimum norm
    #v[k+1]./=exp(1im*angle(v[k+1][1]))
    #res=L(0,0)*v[k+1]+(r+λ[k+1]*L(1,0)*v0)
    #res=LinearAlgebra.norm(res)

    #normalization
    c=0.0+0.0im
    println("order: $k, time: $(time()-zeit)")
    for l=1:k-1
        c-=.5*v[l+1]'*v[k-l+1]
        #c+=v[l+1]'*v[k-l+1]

    end
    v[k+1]+=c*v[1]
    #v[k+1]-=v[1]*(c/2.0+v[1]'*v[k+1])



  end
  return λ,v
end

# function perturb_fast(L,N,v0,v0Adj)
#   #normalize
#   v0/=sqrt(v0'*v0)
#   v0Adj/=v0Adj'*L(1,0)*v0
#   λ=Array{ComplexF64}(undef,N+1)
#   v=Array{Array{Complex{Float64}},1}(undef,N+1)
#   v[1]=v0
#   dim=size(v0)[1]
#   L00=L(0,0)
#   sparse = SparseArrays.issparse(L00)
#   L00=SparseArrays.lu(L(0,0),check=false) #L(0,0) is singular but in most of the cases were are lucky and a LU factorization exists
#   if sparse && L00.status!=0 #lu failed... TODO: implement LU check for dense arrays
#      L00=SparseArrays.qr(L(0,0)) #we try a QR factorization
#   end
#   #TODO consider using qr for all cases
#
#
#   for k=1:N
#     r=zeros(ComplexF64,dim)
#
#     for (m,n) in MN[k+1] #MN is ein dictionary mit schlüsseln (m,n)
#       w=zeros(ComplexF64,dim)
#       for mu in MN[k+1][(m,n)]
#         coeff=1
#         for (g,mu_g) in enumerate(mu)
#           coeff*=λ[g+1]^mu_g #+1 because indexing starts with zero
#           w+=v[k-n-weigh(mu)+1]*multinomcoeff(mu)*coeff
#         end
#       end
#       r+=L(m,n)*w
#     end
#     #v[k+1]= IterativeSolvers.gmres(L(0,0),-(r+λ[k+1]*L(1,0)*v0)) This is crap!!!
#
#     λ[k+1]=-v0Adj'*r/(v0Adj'*L(1,0)*v0)
#     #v[k+1]=L(0,0)\-(r+λ[k+1]*L(1,0)*v0) #This works but LU factorization should be done outside of the loop
#     v[k+1]=L00\-(r+λ[k+1]*L(1,0)*v0) #yeah !
#     res=L(0,0)*v[k+1]+(r+λ[k+1]*L(1,0)*v0)
#   end
#   return λ,v
# end

function perturb_norm(L,N,v0,v0Adj)
  #normalize
  zeit=time()
  Y=-L.terms[end].coeff #TODO: check for __aux__
  Y=convert(SparseArrays.SparseMatrixCSC{ComplexF64,Int},Y)
  v0/=sqrt(v0'*Y*v0)
  v0Adj=SparseArrays.lu(Y)\v0Adj #TODO: this step is highly redundant
  v0Adj/=v0Adj'*Y*L(1,0)*v0
  λ=Array{ComplexF64}(undef,N+1)
  v=Array{Array{Complex{Float64}},1}(undef,N+1)
  v[1]=v0
  dim=size(v0)[1]
  L00=L(0,0)
  sparse = SparseArrays.issparse(L00)
  L00=SparseArrays.lu(L(0,0),check=false) #L(0,0) is singular but in most of the cases were are lucky and a LU factorization exists
  if sparse && L00.status!=0 #lu failed... TODO: implement LU check for dense arrays
     L00=SparseArrays.qr(L(0,0)) #we try a QR factorization
  end
  #TODO consider using qr for all cases
  pack=(@__DIR__)*"/compressed_perturbation_data/"
  for k=1:N
    dir="$k/"
    r=zeros(ComplexF64,dim)
    for m=0:k
      for n=0:k-m
        if m==n==0 || k==m==1
          continue
        end
        fname="$(m)_$(n)"
        w=zeros(ComplexF64,dim)
        open(pack*dir*fname,"r") do file
          for str in eachline(file)
            #str=readline(file)
            mu=part2mult(parse.(Int64,split(str)))

            coeff=1
            for (g,mu_g) in enumerate(mu)
              coeff*=λ[g+1]^mu_g #+1 because indexing starts with zero
            end
            w+=v[k-n-weigh(mu)+1]*multinomcoeff(mu)*coeff
            end
        end
        r+=L(m,n)*w
      end
    end
    #println("##########")
    #println("r: $r")
    #v[k+1]= IterativeSolvers.gmres(L(0,0),-(r+λ[k+1]*L(1,0)*v0)) This is crap!!!
    λ[k+1]=-v0Adj'*Y*r/(v0Adj'*Y*L(1,0)*v0)
    #println("λ: $(λ[k+1])")
    #v[k+1]=L(0,0)\-(r+λ[k+1]*L(1,0)*v0) #This works but LU factorization should be done outside of the loop
    #println("rhs: $(r+λ[k+1]*L(1,0)*v0)")
    v[k+1]=L00\-(r+λ[k+1]*L(1,0)*v0) #yeah !
    v[k+1]-=(v0'*Y*v[k+1])*v0 #if v[k+1] is orthogonal to v0 it features minimum norm
    #v[k+1]./=exp(1im*angle(v[k+1][1]))
    #res=L(0,0)*v[k+1]+(r+λ[k+1]*L(1,0)*v0)
    #res=LinearAlgebra.norm(res)

    #normalization
    c=0.0+0.0im
    println("order: $k, time: $(time()-zeit)")
    for l=1:k-1
        c-=.5*v[l+1]'*Y*v[k-l+1]
        #c+=v[l+1]'*v[k-l+1]

    end
    v[k+1]+=c*v[1]
    #v[k+1]-=v[1]*(c/2.0+v[1]'*v[k+1])



  end
  return λ,v
end
