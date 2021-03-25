function picard(L,z;maxiter=10,tol=0.,relax=1.,n_eig_val=3,v0=[],output=true)
    if output
        println("Launching Picard...")
        println("Iter    Res:     dz:     z:")
        println("----------------------------------")
    end
    z0=complex(Inf)
    n=0
    flag=itsol_converged
    if v0==[]
        v0=ones(ComplexF64,size(L(0))[1])
    end
    M=L(1,oplist=["M"],in_or_ex=true) #get mass matrix
    try
        while abs(z-z0)>tol && n< maxiter
            if output; println(n,"\t\t",abs(z-z0),"\t", z );flush(stdout); end
            z0=z
            X=L(z0,oplist=["M","__aux__"]) #TODO put aux operator as default to linopfamily call and definition
            z,v0=Arpack.eigs(-X,M,sigma=z0,v0=v0,nev=n_eig_val)
            z=sqrt.(z)
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
            println("...aborted Picard!")
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
            println("...finished Picard!")
            println("#####################")
            println(" Picard results ")
            println("#####################")
            println("Number of steps: ",n)
            println("Last step parameter variation:",abs(z0-z))
            println("Eigenvalue:",z)
            println("Eigenvalue/(2*pi):",z/2/pi)
        end
    end
    return Solution(L.params,v0,[],L.eigval), n, flag
end
