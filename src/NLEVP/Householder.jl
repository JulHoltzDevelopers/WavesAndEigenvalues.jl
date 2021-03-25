import Arpack,LinearAlgebra#, NLsolve
#TODO: Degeneracy
# function householder_update(f,a=1)
#     order=length(f)-1
#     a=1.0/a
#
#     if order==1   #aka Newton's method
#         dz=-f[1]/(a*f[2])
#     elseif order==2 #aka Halley's method
#         dz=-2*f[1]*f[2]/((a+1)*f[2]^2-f[1]*f[3])
#     elseif order==3
#         dz=-3*f[1]*((a+1)*f[2]^2-f[1]*f[3]) / ((a^2 + 3*a + 2) *f[2]^3 - 3* (a + 1)* f[1]* f[2] *f[3] + f[1]^2* f[4])
#     elseif order==4
#         dz=-4*f[1]*((a^2 + 3*a + 2) *f[2]^3 - 3* (a + 1)* f[1]* f[2] *f[3] + f[1]^2* f[4]) / (-6* (a^2 + 3 *a + 2) *f[1]* f[2]^2* f[3] + (a^3 + 6 *a^2 + 11* a + 6)* f[2]^4 + f[1]^2 *(3* (a + 1)*f[3]^2 - f[1]* f[5]) + 4* (a + 1)*f[1]^2* f[4] *f[2])
#     else
#         dz=-5*f[1]*(-6* (a^2 + 3 *a + 2) *f[1]* f[2]^2* f[3] + (a^3 + 6 *a^2 + 11* a + 6)* f[2]^4 + f[1]^2 *(3* (a + 1)*f[3]^2 - f[1]* f[5]) + 4* (a + 1)*f[1]^2* f[4] *f[2]) / ((a + 1)*(a + 2)*(a + 3)*(a + 4)* f[2]^5 + 10* (a + 1)*(a + 2)*f[1]^2*f[4]*f[2]^2-10*(a + 1)*(a + 2)*(a + 3)*f[1]*f[2]^3*f[3]+f[1]^3*(f[1]*f[6] - 10 *(a + 1) *f[4]* f[3]) + 5 *(a + 1) *f[1]^2*f[2]* (3* (a + 2)* f[3]^2 - f[1]*f[5]))
#     end
#     return dz
# end

function householder_update(f)
    order=length(f)-1
    if order==1
        dz=-f[1]/f[2]
    elseif order==2
        dz=-f[1]*f[2]/ (f[2]^2-0.5*f[1]*f[3])
    elseif order==3
        dz=-(6*f[1]*f[2]^2-3*f[1]^2*f[3]) / (6*f[2]^3-6*f[1]*f[2]*f[3]+f[1]^2*f[4])
    elseif order == 4
        dz=-(4 *f[1] *(6* f[2]^3 - 6 *f[1] *f[2]* f[3] + f[1]^2* f[4]))/(24 *f[2]^4 - 36* f[1] *f[2]^2 *f[3] + 6*f[1]^2* f[3]^2 + 8*f[1]^2* f[2]* f[4] - f[1]^3* f[5])
    else
        dz=(5 *f[1] *(24 *f[2]^4 - 36* f[1] *f[2]^2* f[3] + 6 *f[1]^2* f[3]^2 + 8* f[1]^2* f[2]* f[4] - f[1]^3* f[5]))/(-120* f[2]^5 + 240* f[1]* f[2]^3* f[3] - 60* f[1]^2* f[2]^2* f[4] + 10* f[1]^2* f[2]* (-9* f[3]^2 + f[1]* f[5]) + f[1]^3* (20* f[3]* f[4] - f[1]* f[6]))
    end
    return dz
end
"""
    sol::Solution, n, flag = householder(L::LinearOperatorFamily, z; <keyword arguments>)

Use a Householder method to iteratively find an eigenpar of `L`, starting the the iteration from `z`.

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
- `output::Bool`: Toggle printing online information.

# Returns
- `sol::Solution`
- `n::Integer`: Number of perforemed iterations
- `flag::Integer`: flag reporting the success of the method. Most important values are `1`: method converged, `0`: convergence might be slow, `-1`:maximum number of iteration has been reached. For other error codes see the source code.

# Notes
Housholder methods are a generalization of Newton's method. If order=1 the Housholder method is identical to Newton's method. The solver then reduces to the "generalized Rayleigh Quotient iteration" presented in [1]. If no relaxation is used (`relax == 1`), the convergence rate is of `order+1`. With relaxation (`relax != 1`) the convergence rate is 1.
Thus, with a higher order less iterations will be necessary. However, the computational time must not necessarily improve nor do the convergence properties. Anyway, if the method converges, the error in the eigenvalue is bounded above by `tol`. For more details on the solver, see the thesis [2].

# References
[1] P. Lancaster, A Generalised Rayleigh Quotient Iteration for Lambda-Matrices,Arch. Rational Mech Anal., 1961, 8, p. 309-322, https://doi.org/10.1007/BF00277446

[2] G.A. Mensah, Efficient Computation of Thermoacoustic Modes, Ph.D. Thesis, TU Berlin, 2019

See also: [`beyn`](@ref)
"""
function householder(L,z;maxiter=10,tol=0.,relax=1.,lam_tol=Inf,order=1,nev=1,v0=[],v0_adj=[],output=true)
    if output
        println("Launching Householder...")
        println("Iter    Res:     dz:     z:")
        println("----------------------------------")
    end
    z0=complex(Inf)
    lam=float(Inf)

    n=0
    active=L.active #IDEA: better integrate activity into householder
    mode=L.mode #IDEA: same for mode
    if v0==[]
        v0=ones(ComplexF64,size(L(0))[1])
    end
    if v0_adj==[]
        v0_adj=conj.(v0)#ones(ComplexF64,size(L(0))[1])
    end



    flag=1
    M=-L.terms[end].coeff
    try
        while  abs(z-z0)>tol   && n<maxiter #&& abs(lam)>lam_tol
            if output; println(n,"\t\t",abs(lam),"\t",abs(z-z0),"\t", z );flush(stdout); end
            z0=z
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
            L.active=[L.auxval,L.eigval]
            for i in 1:nev
                L.params[L.auxval]=lam[i]
                sol=Solution(L.params,v[:,i],v_adj[:,i],L.auxval)
                perturb!(sol,L,L.eigval,order,mode=:householder)
                f=[factorial(idx-1)*coeff for (idx,coeff) in enumerate(sol.eigval_pert[Symbol("$(string(L.eigval))/Taylor")])]
                dz=householder_update(f) #TODO: implement multiplicity
                #print(z+dz)
                push!(delta_z,dz)
            end
            L.active=[L.eigval]
            indexing=sortperm(delta_z, by=abs)
            lam=lam[indexing[1]]
            L.params[L.auxval]=lam #TODO remove this from the loop body
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
            println("...aborted Householder!")
        end

        flag=-2
        if typeof(excp) <: Arpack.ARPACKException
           flag=-4
           if excp==Arpack.ARPACKException(-9999)
               flag=-9999
           end
       elseif excp== LinearAlgebra.SingularException(0) #This means that the solution is so good that L(z) cannot be LU factorized...TODO: implement other strategy into perturb
           flag=-6
           L.params[L.eigval]=z
        end
    end
    if flag==1
        L.params[L.eigval]=z
        if output;  println(n,"\t\t",abs(lam),"\t",abs(z-z0),"\t", z ); end


        if n>=maxiter
            flag=-1
            if output; println("Warning: Maximum number of iterations has been reached!");end

        elseif abs(lam)<=lam_tol
            flag=1
            if output; println("Solution has converged!"); end
        elseif abs(z-z0)<=tol
            flag=0
            if output; println("Warning: Slow convergence!"); end
        elseif isnan(z)
            flag=-5
            if output; println("Warning: computer arithmetics problem. Eigenvalue is NaN"); end
        else
            if output; println("Warning: This should not be possible....\n If you can read this contact GAM!");end
            flag=-3
            println(z)
        end

        if output
            println("...finished Householder!")
            println("#####################")
            println(" Householder results ")
            println("#####################")
            println("Number of steps: ",n)
            println("Last step parameter variation:",abs(z0-z))
            println("Auxiliary eigenvalue λ residual (rhs):", abs(lam))
            println("Eigenvalue:",z)
            println("Eigenvalue/(2*pi):",z/2/pi)
        end
    end
    L.active=active
    L.mode=mode
    #normalization
    v0/=sqrt(v0'*M*v0)
    v0_adj/=conj(v0_adj'*L(L.params[L.eigval],1)*v0)
    return Solution(L.params,v0,v0_adj,L.eigval), n, flag
end

##PAde solver
function poly_roots(p)
  N=length(p)-1
  C=zeros(ComplexF64,N,N)
  for i=2:N
    C[i,i-1]=1
  end
  C[:,N].=-p[1:N]./p[N+1]
  return LinearAlgebra.eigvals(C)
end

function padesolve(L,z;maxiter=10,tol=0.,relax=1.,lam_tol=Inf,order=1,nev=1,v0=[],v0_adj=[],output=true,num_order=1)
    if output
        println("Launching Pade solver...")
        println("Iter    Res:     dz:     z:")
        println("----------------------------------")
    end
    z0=complex(Inf)
    lam=float(Inf)
    lam0=float(Inf)

    n=0
    active=L.active #IDEA: better integrate activity into householder
    mode=L.mode #IDEA: same for mode
    if v0==[]
        v0=ones(ComplexF64,size(L(0))[1])
    end
    if v0_adj==[]
        v0_adj=conj.(v0)#ones(ComplexF64,size(L(0))[1])
    end

    flag=1
    M=-L.terms[end].coeff
    try
        while  abs(z-z0)>tol   && n<maxiter #&& abs(lam)>lam_tol
            if output; println(n,"\t\t",abs(lam),"\t",abs(z-z0),"\t", z );flush(stdout); end

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
            println("...aborted Householder!")
        end

        flag=-2
        if typeof(excp) <: Arpack.ARPACKException
           flag=-4
           if excp==Arpack.ARPACKException(-9999)
               flag=-9999
           end
       elseif excp== LinearAlgebra.SingularException(0) #This means that the solution is so good that L(z) cannot be LU factorized...TODO: implement other strategy into perturb
           flag=-6
           L.params[L.eigval]=z
        end
    end
    if flag==1
        L.params[L.eigval]=z
        if output;  println(n,"\t\t",abs(lam),"\t",abs(z-z0),"\t", z ); end


        if n>=maxiter
            flag=-1
            if output; println("Warning: Maximum number of iterations has been reached!");end

        elseif abs(lam)<=lam_tol
            flag=1
            if output; println("Solution has converged!"); end
        elseif abs(z-z0)<=tol
            flag=0
            if output; println("Warning: Slow convergence!"); end
        elseif isnan(z)
            flag=-5
            if output; println("Warning: computer arithmetics problem. Eigenvalue is NaN"); end
        else
            if output; println("Warning: This should not be possible....\n If you can read this contact GAM!");end
            flag=-3
            println(z)
        end

        if output
            println("...finished Householder!")
            println("#####################")
            println(" Householder results ")
            println("#####################")
            println("Number of steps: ",n)
            println("Last step parameter variation:",abs(z0-z))
            println("Auxiliary eigenvalue λ residual (rhs):", abs(lam))
            println("Eigenvalue:",z)
            println("Eigenvalue/(2*pi):",z/2/pi)
        end
    end
    L.active=active
    L.mode=mode
    #normalization
    v0/=sqrt(v0'*M*v0)
    v0_adj/=conj(v0_adj'*L(L.params[L.eigval],1)*v0)
    return Solution(L.params,v0,v0_adj,L.eigval), n, flag
end
