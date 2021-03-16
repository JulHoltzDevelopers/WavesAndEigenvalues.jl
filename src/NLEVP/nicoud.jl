function nicoud(L,z;maxiter=10,tol=0.,relax=1.,n_eig_val=3,v0=[],output=true)
    if output
        println("Launching Nicoud...")
        println("Iter    Res:     dz:     z:")
        println("----------------------------------")
    end
    z0=complex(Inf)
    n=0
    flag=itsol_converged

    M=L(1,oplist=["M"],in_or_ex=true) #get mass matrix
    K=L(1,oplist=["K"],in_or_ex=true)
    C=L(1,oplist=["C"],in_or_ex=true)
    d=size(M)[1] #TODO: implement a dimension method for LinearOperatorFamily
    I=SparseArrays.sparse(SparseArrays.I,d,d)
    O=SparseArrays.spzeros(d,d)
    Y=[I O; O M]
    if v0==[]
        v0=ones(ComplexF64,d)
    end
    v0=[v0; z*v0]
    try
        while abs(z-z0)>tol && n< maxiter
            if output; println(n,"\t\t",abs(z-z0),"\t", z );flush(stdout); end
            z0=z
            Q=L(z,oplist=["Q"],in_or_ex=true)
            X=[O -I; K+Q C]
            z,v0=Arpack.eigs(-X,Y,sigma=z0,v0=v0,nev=n_eig_val)
            indexing=sortperm(z.-z0,by=abs)
            z,v0=z[indexing[1]],v0[:,indexing[1]]
            z=z0+relax*(z-z0)
            #TODO consider relaxation in v0
            n+=1
        end
    catch excp
        if output
            println("Error occured:")
            println(excp)
            println(typeof(excp))
            println("...aborted Nicoud!")
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
        if output;  println(n,"\t\t",abs(z-z0),"\t", z ); end


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
        end
        if output
            println("...finished Nicoud!")
            println("#####################")
            println(" Nicoud results ")
            println("#####################")
            println("Number of steps: ",n)
            println("Last step parameter variation:",abs(z0-z))
            println("Eigenvalue:",z)
            println("Eigenvalue/(2*pi):",z/2/pi)
        end
    end
    return Solution(L.params,v0[1:d],[],L.eigval), n, flag
end
