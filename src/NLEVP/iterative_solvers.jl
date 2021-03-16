
## status flags for iterative solvers
#found a solution  = zero
const itsol_converged = 0
# warnings (positive numbers)
const itsol_maxiter = 1
const itsol_slow_convergence = 2
#errors (negative numbers)
const itsol_impossible = -1
const itsol_singular_exception = -2
const itsol_arpack_exception = -3
const itsol_isnan = -4
const itsol_unknown = -5
const itsol_arpack_9999 = -9999
"""
    msg=decode_error_flag(flag::Int)

Decode error flag `flag` from an iterative solver into a string `msg`.

See also: [`inveriter`](@ref), [`lancaster`](@ref), [`mslp`](@ref), [`rf2s`](@ref), [`traceiter`](@ref)
"""
function decode_error_flag(flag::Int)
    if flag == itsol_converged
        msg = "Solution converged, everythink OK!"
    elseif flag == itsol_maxiter
        msg = "Warning: Maximum number of iterations has been reached!"
    elseif flag == itsol_slow_convergence
        msg = "Warning: Slow progress!"
    elseif flag == itsol_impossible
        msg == "Error: This error should be impossible. Please, contact the package developers!"
    elseif flag == itsol_singular_exception.
        msg == "Error: Singular Exception!"
    elseif flag == itsol_arpack_exception
        msg == "Error: Arpack exception!"
    elseif flag == itsol_arpack_999
        msg == "Error: Arpack -9999 error!"
    elseif flag itsol_unknown
        msg == "Error: Unknown error ocurred!"
    else
        msg == "Unknown flag code."
    end

    return msg
end



## method of successive linear problems
"""
    sol,n,flag = mslp(L,z;maxiter=10,tol=0.,relax=1.,lam_tol=Inf,order=1,nev=1,v0=[],v0_adj=[],num_order=1,scale=1.0+0.0im,output=false)

Deploy the method of successive linear problems for finding an eigentriple of `L`.

# Arguments
- `L::LinearOperatorFamily`: Definition of the nonlinear eigenvalue problem.
- `z`: Initial guess for the eigenvalue.
- `maxiter::Integer=10`: Maximum number of iterations.
- `tol=0`: Absolute tolerance to trigger the stopping of the iteration. If the difference of two consecutive iterates is `abs(z0-z1)<tol` the iteration is aborted.
- `relax=1`: relaxation parameter
- `lam_tol=Inf`: tolerance for the auxiliary eigenvalue to test convergence. The default is infinity, so there is effectively no test.
- `order::Integer=1`: Order of the Householder method, Maximum is 5
- `nev::Integer=1`: Number of Eigenvalues to be searched for in intermediate ARPACK calls.
- `v0::Vector`: Initial vector for Krylov subspace generation in ARPACK calls. If not provided the vector is initialized with ones.
- `v0_adj::Vector`: Initial vector for Krylov subspace generation in ARPACK calls. If not provided the vector is initialized with `v0`.
- `num_order::Int=1`: order of the numerator polynomial when forming the Padé approximant.
- `scale::Number=1`: scaling factor applied to the tolerance (and the eigenvalue output when `output==true`).
- `output::Bool=false`: Toggle printing online information.

# Returns
- `sol::Solution`
- `n::Integer`: Number of perforemed iterations
- `flag::Integer`: flag reporting the success of the method. Most important values are `1`: method converged, `0`: convergence might be slow, `-1`:maximum number of iteration has been reached. For other error codes see the source code.

# Notes
This is a variation of the method of succesive linear Problems [1] as published in Algorithm 3 in [2].
It treats the eigenvalue `ω` as a parameter in a linear eigenvalue problem `L(ω)p=λYp` and computes the root of the implicit relation `λ(ω)==0`.

The default values "order==1" and "num_order==1" will result in a Newton-iteration for finding the roots. Choosing a higher order will result in
Housholder methods as a generalization of Newton's method. If no relaxation is used (`relax == 1`), the convergence rate is of `order+1`. With relaxation (`relax != 1`) the convergence rate is 1.
Thus, with a higher order less iterations will be necessary. However, the computational time must not necessarily improve nor do the convergence properties. Anyway, if the method converges, the error in the eigenvalue is bounded above by `tol`. For more details on the solver, see the thesis [2].

Householder methods are based on expanding the relation `λ=λ(ω)` into a `[1/order-1]`-Padé approximant. The numerator order may be increased using the keyword "num_order". This feature is experimental.

# References
[1] A.Ruhe, Algorithms for the non-linear eigenvalue problem, SIAM J. Numer. Anal. ,10:674–689, 1973.

[2] V. Mehrmann and H. Voss (2004), Nonlinear eigenvalue problems: A challenge for modern eigenvalue methods, GAMM-Mitt. 27, 121–152.

[3] G.A. Mensah, Efficient Computation of Thermoacoustic Modes, Ph.D. Thesis, TU Berlin, 2019 [doi:10.14279/depositonce-8952 ](http://dx.doi/10.14279/depositonce-8952)

See also: [`beyn`](@ref), [`inveriter`](@ref), [`lancaster`](@ref), [`rf2s`](@ref), [`traceiter`](@ref), [`decode_error_flag`](@ref)
"""
function mslp(L,z;maxiter=10,tol=0.,relax=1.,lam_tol=Inf,order=1,nev=1,v0=[],v0_adj=[],num_order=1, scale=1, output=true)
    if output
        println("Launching MSLP solver...")
        if scale!=1
            println("scale: $scale")
        end
        println("Iter   dz:     z:")
        println("----------------------------------")
    end
    z0=complex(Inf)
    lam=float(Inf)
    lam0=float(Inf)
    z*=scale
    tol*=scale

    n=0
    active=L.active #IDEA: better integrate activity into householder
    mode=L.mode #IDEA: same for mode
    if v0==[]
        v0=ones(ComplexF64,size(L(0))[1])
    end
    if v0_adj==[]
        v0_adj=conj.(v0)#ones(ComplexF64,size(L(0))[1])
    end

    flag=itsol_converged
    if L.terms[end].operator != "__aux__"
        push!(L,Term(-LinearAlgebra.I,(pow1,),((:__aux__,),),"__aux__","__aux__"))
        L.auxval=:__aux__
        #IDEA: possibly remove term at the end of the function
    end
    M=-L.terms[end].coeff
    try
        while  abs(z-z0)>tol   && n<maxiter #&& abs(lam)>lam_tol
            if output; println(n,"\t\t","\t",abs(z-z0)/scale,"\t", z/scale );flush(stdout); end

            L.params[L.eigval]=z
            L.params[L.auxval]=0
            A=L(z)
            lam,v = Arpack.eigs(A,M,nev=nev, sigma = 0,v0=v0)
            lam_adj,v_adj = Arpack.eigs(A',M', nev=nev, sigma = 0,v0=v0_adj)
            #TODO: consider constdouton
            #TODO: multiple eigenvalues
            indexing=sortperm(lam, by=abs)
            lam=lam[indexing]
            v=v[:,indexing]
            indexing=sortperm(lam_adj, by=abs)
            lam_adj=lam_adj[indexing]
            v_adj=v_adj[:,indexing]
            delta_z =[]
            back_delta_z=[]
            L.active=[L.auxval,L.eigval]
            #println("#############")
            #println("lam0:$lam0 ")
            for i in 1:nev
                L.params[L.auxval]=lam[i]
                sol=Solution(L.params,v[:,i],v_adj[:,i],L.auxval)
                perturb!(sol,L,L.eigval,order,mode=:householder)
                coeffs=sol.eigval_pert[Symbol("$(L.eigval)/Taylor")]
                num,den=pade(coeffs,num_order,order-num_order)

                #forward calculation (has issues with multi-valuedness)
                roots=poly_roots(num)
                indexing=sortperm(roots, by=abs)
                dz=roots[indexing[1]]
                push!(delta_z,dz)
                #poles=sort(poly_roots(den),by=abs)
                #poles=poles[1]
                #println(">>>$i<<<")
                #println("$coeffs")
                #println("residue:$(LinearAlgebra.norm((A-lam[i]*M)*v[:,i])) and $(LinearAlgebra.norm((A'-lam_adj[i]*M')*v_adj[:,i]))")
                #println("poles: $(poles+z) r:$(abs(poles))")
                #backward check(for solving multi-valuedness problems by continuity)
                if z0!=Inf
                    back_lam=polyval(num,z0-z)/polyval(den,z0-z)
                    #println("back:$back_lam")
                    #println("root: $(z+dz)")
                    #estm_lam=polyval(num,dz)/polyval(den,dz)
                    #println("estm. lam: $estm_lam")
                    back_lam=lam0-back_lam
                    push!(back_delta_z,back_lam)
                end
            end
            L.active=[L.eigval]
            if z0!=Inf
                indexing=sortperm(back_delta_z, by=abs)
            else
                indexing=sortperm(delta_z, by=abs)
            end
            lam=lam[indexing[1]]
            L.params[L.auxval]=lam #TODO remove this from the loop body
            z0=z
            lam0=lam
            z=z+relax*delta_z[indexing[1]]
            v0=(1-relax)*v0+relax*v[:,indexing[1]]
            v0_adj=(1-relax)*v0_adj+relax*v_adj[:,indexing[1]]
            n+=1
        end

    catch excp
        if output
            println("Error occured:")
            println(excp)
            println(typeof(excp))
            println("...aborted MSLP!")
        end

        flag=itsol_unknown
        if typeof(excp) <: Arpack.ARPACKException
           flag=itsol_arpack_exception
           if excp==Arpack.ARPACKException(-9999)
               flag=itsol_arpack_9999
           end
       elseif excp== LinearAlgebra.SingularException(0) #This means that the solution is so good that L(z) cannot be LU factorized...TODO: implement other strategy into perturb
           flag=itsol_singular_exception
           L.params[L.eigval]=z
        end
    end
    if flag==itsol_converged
        L.params[L.eigval]=z
        if output;  println(n,"\t\t",abs(lam),"\t",abs(z-z0),"\t", z ); end


        if n>=maxiter
            flag=itsol_maxiter
            if output; println("Warning: Maximum number of iterations has been reached!");end

        elseif abs(lam)<=lam_tol
            flag=itsol_converged
            if output; println("Solution has converged!"); end
        elseif abs(z-z0)<=tol
            flag=itsol_slow_convergence
            if output; println("Warning: Slow convergence!"); end
        elseif isnan(z)
            flag=itsol_isnan
            if output; println("Warning: computer arithmetics problem. Eigenvalue is NaN"); end
        else
            if output; println("Warning: This should not be possible....\n If you can read this contact GAM!");end
            flag=itsol_impossible
            println(z)
        end

        if output
            println("...finished MSLP!")
            println("#####################")
            println(" Results ")
            println("#####################")
            println("Number of steps: ",n)
            println("Last step parameter variation:",abs(z0-z))
            println("Auxiliary eigenvalue $(L.auxval) residual (rhs):", abs(lam))
            println("Eigenvalue:",z/scale)
        end
    end
    L.active=active
    L.mode=mode
    #normalization
    v0/=sqrt(v0'*M*v0)
    v0_adj/=conj(v0_adj'*L(L.params[L.eigval],1)*v0)
    return Solution(L.params,v0,v0_adj,L.eigval), n, flag
end

## inverse iterations
"""
    sol,n,flag = inveriter(L,z;maxiter=10,tol=0., relax=1., x0=[],v=[], output=true)

Deploy inverse iteration for finding an eigenpair of `L`.

# Arguments
- `L::LinearOperatorFamily`: Definition of the nonlinear eigenvalue problem.
- `z`: Initial guess for the eigenvalue.
- `maxiter::Integer=10`: Maximum number of iterations.
- `tol=0`: Absolute tolerance to trigger the stopping of the iteration. If the difference of two consecutive iterates is `abs(z0-z1)<tol` the iteration is aborted.
- `relax=1`: relaxation parameter
- `x0::Vector`: initial guess for eigenvector. If not provided it is all ones.
- `v0::Vector`: normalization vector. If not provided it is all ones.
- `output::Bool`: Toggle printing online information.

# Returns
- `sol::Solution`
- `n::Integer`: Number of perforemed iterations
- `flag::Integer`: flag reporting the success of the method. Most important values are `1`: method converged, `0`: convergence might be slow, `-1`:maximum number of iteration has been reached. For other error codes see the source code.

# Notes
The algorithm uses Newton-Raphson itreations on the vector level for iteratively solving the problem
`L(λ)x == 0 ` and `v'x == 0`.
The implementation is based on Algorithm 1 in [1]

# References
[1] V. Mehrmann and H. Voss (2004), Nonlinear eigenvalue problems: A challenge for modern eigenvalue methods, GAMM-Mitt. 27, 121–152.

See also: [`beyn`](@ref), [`lancaster`](@ref), [`mslp`](@ref), [`rf2s`](@ref), [`traceiter`](@ref), [`decode_error_flag`](@ref)
"""
function inveriter(L,z;maxiter=10,tol=0., relax=1., x0=[],v=[], output=true)
  if output
      println("Launching inverse iteration...")
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
  flag=itsol_converged
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
        println("...aborted inverse iteration!")
    end

    flag=itsol_unknown
  end

  #convergence checks
  if flag==itsol_converged
    if output;  println(n,"\t\t",abs(z-z0),"\t", z );flush(stdout); end
    if n>=maxiter
      flag=itsol_maxiter
      if output; println("Warning: Maximum number of iterations has been reached!");end
    elseif abs(z-z0)<=tol
      flag=itsol_converged
      if output; println("Solution has converged!"); end
    else
      if output; println("Warning: This should not be possible....\n If you can read this contact GAM!");end
      flag=itsol_impossible
      println(z)
    end
  end

  y=[]

  return  Solution(L.params,x0,y,L.eigval,L.auxval), n , flag
end

## Lancaster's Rayleigh-quotient iteration
"""
    sol,n,flag = lancaster(L,z;maxiter=10,tol=0., relax=1., x0=[],y0=[], output=true)

Deploy Lancaster's Rayleigh-quotient iteration for finding an eigenvalue of `L`.

# Arguments
- `L::LinearOperatorFamily`: Definition of the nonlinear eigenvalue problem.
- `z`: Initial guess for the eigenvalue.
- `maxiter::Integer=10`: Maximum number of iterations.
- `tol=0`: Absolute tolerance to trigger the stopping of the iteration. If the difference of two consecutive iterates is `abs(z0-z1)<tol` the iteration is aborted.
- `relax=1`: relaxation parameter
- `x0::Vector`: initial guess for right iteration vector. If not provided it is all ones.
- `y0::Vector`: initial guess for left iteration vector. If not provided it is all ones.
- `output::Bool`: Toggle printing online information.

# Returns
- `sol::Solution`
- `n::Integer`: Number of perforemed iterations
- `flag::Integer`: flag reporting the success of the method. Most important values are `1`: method converged, `0`: convergence might be slow, `-1`:maximum number of iteration has been reached. For other error codes see the source code.

# Notes
The algorithm is a generalization of Rayleigh-quotient-iteration by Lancaster [1].

# References
[1] P. Lancaster, A Generalised Rayleigh Quotient Iteration for Lambda-Matrices,Arch. Rational Mech Anal., 1961, 8, p. 309-322, https://doi.org/10.1007/BF00277446

See also: [`beyn`](@ref), [`inveriter`](@ref), [`mslp`](@ref), [`rf2s`](@ref), [`traceiter`](@ref), [`decode_error_flag`](@ref)
"""
function lancaster(L,z;maxiter=10,tol=0., relax=1., x0=[],y0=[], output=true)
  if output
      println("Launching Lancaster's Rayleigh-quotient iteration...")
  end
  if x0==[]
    x0=ones(ComplexF64,size(L(0))[1])
  end
  if y0==[]
    y0=ones(ComplexF64,size(L(0))[1])
  end

  z0=complex(Inf)
  n=0
  flag=itsol_converged
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
        println("...aborted Lancaster's method!")
    end

    flag=itsol_unknown
  end

  #convergence checks
  if flag==itsol_converged
    if output;  println(n,"\t\t",abs(z-z0),"\t", z );flush(stdout); end
    if n>=maxiter
      flag=itsol_maxiter
      if output; println("Warning: Maximum number of iterations has been reached!");end
    elseif abs(z-z0)<=tol
      flag=itsol_converged
      if output; println("Solution has converged!"); end
    else
      if output; println("Warning: This should not be possible....\n If you can read this contact GAM!");end
      flag=itsol_impossible
      println(z)
    end
  end


  return  Solution(L.params,zeros(ComplexF64,length(x0)),[],L.eigval), n , flag
end

## trace iteration
"""
    sol,n,flag = traceiter(L,z;maxiter=10,tol=0.,relax=1.,output=true)

Deploy trace iteration for finding an eigevalue of `L`.

# Arguments
- `L::LinearOperatorFamily`: Definition of the nonlinear eigenvalue problem.
- `z`: Initial guess for the eigenvalue.
- `maxiter::Integer=10`: Maximum number of iterations.
- `tol=0`: Absolute tolerance to trigger the stopping of the iteration. If the difference of two consecutive iterates is `abs(z0-z1)<tol` the iteration is aborted.
- `relax=1`: relaxation parameter
- `output::Bool`: Toggle printing online information.

# Returns
- `sol::Solution`
- `n::Integer`: Number of perforemed iterations
- `flag::Integer`: flag reporting the success of the method. Most important values are `1`: method converged, `0`: convergence might be slow, `-1`:maximum number of iteration has been reached. For other error codes see the source code.

# Notes
The algorithm applies Newton-Raphson iteration for finding the root of the determinant. The updates for the iterates are computed by using trace operations and  Jacobi's formula [1]. The trace operation renders the method slow and inaccurate for medium and large problems.

# References
[1] S. Güttel and F. Tisseur, The Nonlinear Eigenvalue Problem, 2017, http://eprints.ma.man.ac.uk/2531/.

See also: [`beyn`](@ref), [`inveriter`](@ref), [`lancaster`](@ref), [`mslp`](@ref), [`rf2s`](@ref), [`decode_error_flag`](@ref)
"""
function traceiter(L,z;maxiter=10,tol=0.,relax=1.,output=true)
    if output
        println("Launching trace iteration...")
    end
  z0=complex(Inf)
  n=0
  flag=itsol_converged
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
        println("...aborted trace iteration!")
    end
    flag=itsol_unknown
  end


  #convergence checks
  if flag==itsol_converged
    if output;  println(n,"\t\t",abs(z-z0),"\t", z );flush(stdout); end
    if n>=maxiter
      flag=itsol_maxiter
      if output; println("Warning: Maximum number of iterations has been reached!");end
    elseif abs(z-z0)<=tol
      flag=itsol_converged
      if output; println("Solution has converged!"); end
    else
      if output; println("Warning: This should not be possible....\n If you can read this contact GAM!");end
      flag=itsol_impossible
      println(z)
    end
  end

    return Solution(L.params,[],[],L.eigval), n, flag
end

## two-sided Rayleigh-functional iteration
"""
    sol,n,flag = rf2s(L,z;maxiter=10,tol=0., relax=1., x0=[],y0=[], output=true)

Deploy two sided Rayleigh-functional iteration for finding an eigentriple of `L`.

# Arguments
- `L::LinearOperatorFamily`: Definition of the nonlinear eigenvalue problem.
- `z`: Initial guess for the eigenvalue.
- `maxiter::Integer=10`: Maximum number of iterations.
- `tol=0`: Absolute tolerance to trigger the stopping of the iteration. If the difference of two consecutive iterates is `abs(z0-z1)<tol` the iteration is aborted.
- `relax=1`: relaxation parameter
- `x0::Vector`: initial guess for right eigenvector. If not provided it is the first basis vector `[1 0 0 0 ...]`.
- `y0::Vector`: initial guess for left eigenvector. If not provided it is the first basis vector `[1 0 0 0 ...].
- `output::Bool`: Toggle printing online information.

# Returns
- `sol::Solution`
- `n::Integer`: Number of perforemed iterations
- `flag::Integer`: flag reporting the success of the method. Most important values are `1`: method converged, `0`: convergence might be slow, `-1`:maximum number of iteration has been reached. For other error codes see the source code.

# Notes
The algorithm is known to have a cubic convergence rate. Its implementation follows Algorithm 4.9 in [1].

# References
[1] S. Güttel and F. Tisseur, The Nonlinear Eigenvalue Problem, 2017, http://eprints.ma.man.ac.uk/2531/.

See also: [`beyn`](@ref), [`inveriter`](@ref), [`lancaster`](@ref), [`mslp`](@ref), [`traceiter`](@ref), [`decode_error_flag`](@ref)
"""
function rf2s(L,z;maxiter=10,tol=0., relax=1., x0=[],y0=[], output=true)
  if output
      println("Launching two-sided Rayleigh functional iteration...")
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
  flag=itsol_converged
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
        println("...aborted two-sided Rayleigh functional iteration!")
    end
    flag=itsol_unknown
  end
  #convergence checks
  if flag==itsol_converged
    if output;  println(n,"\t\t",abs(z-z0),"\t", z );flush(stdout); end
    if n>=maxiter
      flag=itsol_maxiter
      if output; println("Warning: Maximum number of iterations has been reached!");end
    elseif abs(z-z0)<=tol
      flag=itsol_converged
      if output; println("Solution has converged!"); end
    elseif isnan(z)
        flag=itsol_isnan
        if output; println("Warning: computer arithmetics problem. Eigenvalue is NaN"); end
    else
      if output; println("Warning: This should not be possible....\n If you can read this contact GAM!");end
      flag=itsol_impossible
      println(z)
    end
  end
      return Solution(L.params,x0,y0,L.eigval), n, flag
end
