function mehrmann(L,z;maxiter=10,tol=0., relax=1., x0=[],v=[], output=true)
  if output
      println("Launching Mehrmann...")
  end
  if x0==[]
    x0=ones(ComplexF64,size(L(0))[1])
  end
  if v==[]
    v=ones(ComplexF64,size(L(0))[1])
  end
  #normalize
  x0/=v'*x0

  active=L.active #TODO: better integrate activity into householder

  z0=complex(Inf)
  n=0
  flag=1
  try
    while abs(z-z0)>tol && n<maxiter
      if output;  println(n,"\t\t",abs(z-z0),"\t", z );flush(stdout); end
      z0=z
      u=L(z,0)\(L(z,1)*x0)
      #println(size(v)," ",size(x0)," ",size(u))
      z=z0-(v'*x0)/(v'*u)
      x0=u/(v'*u)
      #println("###",v'x0)
      n+=1
    end
  catch excp
    if output
        println("Error occured:")
        println(excp)
        println(typeof(excp))
        println("...aborted Mehrmann!")
    end

    flag=-2
  end

  #convergence checks
  if flag==1
    if output;  println(n,"\t\t",abs(z-z0),"\t", z );flush(stdout); end
    if n>=maxiter
      flag=-1
      if output; println("Warning: Maximum number of iterations has been reached!");end
    elseif abs(z-z0)<=tol
      flag=1
      if output; println("Solution has converged!"); end
    else
      if output; println("Warning: This should not be possible....\n If you can read this contact GAM!");end
      flag=-3
      println(z)
    end
  end

  # if flag==1
  #   #compute left eigenvector
  #   #TODO: degeneracy
  #   lam,y = Arpack.eigs(L(z)',nev=1, sigma = 0,v0=x0)
  #   #normalize
  #   println("size:::$size(y)")
  #   x0./sqrt(x0'*x0)
  #   y./=conj(y'*L(z,1)*x0)
  # else
  #   y=[]
  # end
  y=[]

  return  Solution(L.params,x0,y,L.eigval), n , flag

end
##
function lancaster(L,z;maxiter=10,tol=0., relax=1., x0=[],y0=[], output=true)
  if output
      println("Launching Lancaster...")
  end
  if x0==[]
    x0=ones(ComplexF64,size(L(0))[1])
  end
  if y0==[]
    y0=ones(ComplexF64,size(L(0))[1])
  end

  z0=complex(Inf)
  n=0
  flag=1
  try
    while abs(z-z0)>tol && n<maxiter
            if output;  println(n,"\t\t",abs(z-z0),"\t", z );flush(stdout); end
      z0=z
      #solve
      ξ=L(z)\x0
      η=L(z)'\y0
      z=z0-(η'*L(z,0)*ξ)/(η'*L(z,1)*ξ)
      n+=1
    end
  catch excp
    if output
        println("Error occured:")
        println(excp)
        println(typeof(excp))
        println("...aborted Lancaster!")
    end

    flag=-2
  end

  #convergence checks
  if flag==1
    if output;  println(n,"\t\t",abs(z-z0),"\t", z );flush(stdout); end
    if n>=maxiter
      flag=-1
      if output; println("Warning: Maximum number of iterations has been reached!");end
    elseif abs(z-z0)<=tol
      flag=1
      if output; println("Solution has converged!"); end
    else
      if output; println("Warning: This should not be possible....\n If you can read this contact GAM!");end
      flag=-3
      println(z)
    end
  end


  return  Solution(L.params,zeros(ComplexF64,length(x0)),[],L.eigval), n , flag
end






#%%
import LinearAlgebra
function juniper(L,z;maxiter=10,tol=0.,relax=1.,output=true)
    if output
        println("Launching Juniper...")
    end
  z0=complex(Inf)
  n=0
  flag=1
  try
    while abs(z-z0)>tol && n<maxiter
      if output;  println(n,"\t\t",abs(z-z0),"\t", z ); end
      z0=z
      #dz=-1/LinearAlgebra.tr(Array(L(z))\Array(L(z,1))) #TODO: better implementation for sparse types
      #QR=SparseArrays.qr(L(z))
      LU=SparseArrays.lu(L(z),check=false)
      L1=L(z,1)
      dz=0.0+0.0im
      for idx=1:size(LU)[1]
        dz+=(LU\Array(L1[:,idx]))[idx]
      end
      dz=-1/dz
      z=z0+relax*dz
      n+=1
    end
  catch excp
    if output
        println("Error occured:")
        println(excp)
        println(typeof(excp))
        println("...aborted Juniper!")
    end
    flag=-2
  end


  #convergence checks
  if flag==1
    if output;  println(n,"\t\t",abs(z-z0),"\t", z );flush(stdout); end
    if n>=maxiter
      flag=-1
      if output; println("Warning: Maximum number of iterations has been reached!");end
    elseif abs(z-z0)<=tol
      flag=1
      if output; println("Solution has converged!"); end
    else
      if output; println("Warning: This should not be possible....\n If you can read this contact GAM!");end
      flag=-3
      println(z)
    end
  end

    return Solution(L.params,[],[],L.eigval), n, flag
end


##

function guettel(L,z;maxiter=10,tol=0., relax=1., x0=[],y0=[], output=true)
  if output
      println("Launching Guettel...")
  end
  if x0==[]
    x0=zeros(ComplexF64,size(L(0))[1])
    x0[1]=1
  end
  if y0==[]
    y0=zeros(ComplexF64,size(L(0))[1])
    y0[1]=1
  end
  #normalize
  x0/=sqrt(x0'*x0)
  y0/=sqrt(y0'*y0)
  z0=complex(Inf)
  n=0
  flag=1
  try
    while abs(z-z0)>tol && n<maxiter
      if output;  println(n,"\t\t",abs(z-z0),"\t", z ); end
      z0=z
      LU=SparseArrays.lu(L(z),check=false)
      x0=LU\(L(z,1)*x0)
      y0=LU'\(L(z,1)'*y0)
      #normalize
      x0/=sqrt(x0'*x0)
      y0/=sqrt(y0'*y0)
      #inner iteration
      idx=0
      z00=complex(Inf)
      while abs(z-z00)>tol && idx<10
        z00=z
        z=z-(y0'*L(z)*x0)/(y0'*L(z,1)*x0)
        idx+=1
      end
      n+=1
    end
  catch excp
    if output
        println("Error occured:")
        println(excp)
        println(typeof(excp))
        println("...aborted Guettel!")
    end
    flag=-2
  end
  #convergence checks
  if flag==1
    if output;  println(n,"\t\t",abs(z-z0),"\t", z );flush(stdout); end
    if n>=maxiter
      flag=-1
      if output; println("Warning: Maximum number of iterations has been reached!");end
    elseif abs(z-z0)<=tol
      flag=1
      if output; println("Solution has converged!"); end
    elseif isnan(z)
        flag=-5
        if output; println("Warning: computer arithmetics problem. Eigenvalue is NaN"); end
    else
      if output; println("Warning: This should not be possible....\n If you can read this contact GAM!");end
      flag=-3
      println(z)
    end
  end
      return Solution(L.params,x0,y0,L.eigval), n, flag
end


## count poles and zeros

function count_eigvals(L,Γ,N;output=false)
  function integrand(z)
    z,w=z
    LU=SparseArrays.lu(L(z),check=false)
    L1=L(z,1)
    dz=0.0+0.0im
    for idx=1:size(LU)[1]
      dz+=(LU\Array(L1[:,idx]))[idx]
    end
    return dz*w
  end

  return gauss(0.0+0.0im,integrand,Γ,N,output)/2/pi/1.0im
end
