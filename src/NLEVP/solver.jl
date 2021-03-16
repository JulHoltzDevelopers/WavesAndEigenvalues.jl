"""
    solve(L,Γ;<keywordarguments>)

Find eigenvalues of linear operator family `L` inside contour `Γ`.

# Arguments
- `L::LinearOperatorFamily`: Linear operator family for which eigenvalues are to be found
- `Γ::List`: List of complex points defining the contour.

# Keywordarguments (optional)
- `Δl::Int=1`: Increment to the search space dimension in Beyn's algorithm
- `N::Int=16`: Number of integration points per edge of `Γ`.
- `tol::Float=1E-8`: Tolerance for local eigenvalue solver.
- `eigvals::Dict=Dict()`: Dictionary of known eigenvalues and corresponding eigenvectors.
- `maxcycles::Int=1`: Number of local correction cycles to global eigenvalue estimate.
- `nev::Int=1`: Number of tested eigenvalues for finding the best local match to the global estimates.
- `max_outer_cycles::Int=1`: Number of outer cycles, i.e., incremental increases to the search space dimension by `Δl`.
- `atol_σ::Float=1E-12`: absolute tolerance for terminating the inner cycle if the maximum singular value is `σ<atol_σ``
- `rtol_σ::Float=1E-8`: relative tolerance for terminating the inner cycle if the the ratio of the initial and the current maximum singular value is `σ/σ0<rtol_σ`
- `loglevel::Int=0`: loglevel value between 0 (no output) and 2 (max output) to control the detail of printed output.

# Returns
- `eigvals::Dict`:List of computed eigenvalues and corresponding eigenvectors.

# Notes
The algorithm is experimental. It startes with a low-dimensional Beyn integral.
From which estimates to the eigenvalues are computed. These estimates are used
to analytically correct Beyn's integral. This step requires no additional
numerical integration and is therefore relatively quick. New local estimates are
computed and if the local solver then found new eigenvalues these are used to
further correct Beyn's integral. In an optional outer loop the entire procedure
can be repeated, after increasing the search space dimension.

See also: [`householder`](@ref), [`beyn`](@ref)
"""
function solve(L::LinearOperatorFamily,Γ;Δl=1,N=16, tol=1E-8,eigvals=Dict(),maxcycles=1,nev=1,max_outer_cycles=1,atol_σ=1E-12,rtol_σ=1E-8,loglevel=0)
    #header
    d=size(L.terms[1].coeff)[1]#system dimension
    A=[]
    l=Δl
    σ0,σ,σmax=0.0,0.0,0.0
    #outer loop
    while l<=max_outer_cycles*Δl
        V=zeros(ComplexF64,d,Δl)
        for ll=1:Δl
            #TODO: randomization
            V[(l-Δl)+ll,ll]=1.0+0.0im
        end
        A=append!(A,[compute_moment_matrices(L,Γ,V,K=1,N=N,output=true),])
        #σmax=max(σmax,maximum(abs.(A[1][:,:,1])))
        if l>Δl #reinitialize σmax in new iteration of outer loop
            Ω,P,Σ=moments2eigs(A,return_σ=true)
            #println(maximum(Σ))
            σmax,σ0,σ=max(σmax,maximum(Σ)),maximum(Σ),0.0
        end
        #update moment matrix with known solutions
        for (ω, val) in eigvals
            w=wn(ω,Γ)
            for ll=1:Δl
                moment=-2*pi*1.0im*w*val[1].v*val[1].v_adj[ll]'
                A[l÷Δl][:,ll,1]+=moment
                A[l÷Δl][:,ll,2]+=ω*moment
            end
        end

        number_of_interior_eigvals=0
        for val in values(eigvals)
            if val[2]
                number_of_interior_eigvals+=1
            end
        end
        neigval=number_of_interior_eigvals
        cycle=0
        while cycle<maxcycles# inner cycle
            cycle+=1
            if loglevel>=1
                flush(stdout)
                println("####cycle$cycle####")
            end
            Ω,P,Σ=moments2eigs(A,return_σ=true)
            #println(maximum(Σ))
            σmax,σ0,σ=max(σmax,σ),σ,maximum(Σ)


            #compute accurate local solutions
            for idx=1:length(Ω)
                ω=Ω[idx]
                if loglevel>=2
                    flush(stdout)
                    println("omega:$ω")
                end
                #remove known vectors from first Krylov vector
                 v0=P[:,idx]
                v0/=sqrt(v0'v0)
                for (key,val) in eigvals
                    v=val[1].v
                    v/=sqrt(v'v)
                    h=v'v0
                    v0-=h*v
                    v0/=sqrt(v0'v0)
                end
                # #TODO: initilaize adjoint vector
                #
                #
                # sol,nn,flag=householder(L,ω,maxiter=10,tol=tol,output=false,order=3,nev=nev,v0=v0)
                sol,nn,flag=mehrmann(L,ω,maxiter=10,tol=tol,output=false,x0=v0,v=v0)
                ω=sol.params[sol.eigval]
                if flag>=0
                    is_new_value=true
                    for val in keys(eigvals)
                        if abs(ω-val)<10*tol
                            is_new_value=false
                            #break
                        end
                    end
                else
                    is_new_value=false
                end
                if loglevel>=2
                    println("conv:$ω flag:$flag new: $is_new_value")
                end
                #update moment matrix with local solutions
                if is_new_value && inpoly(ω,Γ)
                    #R=5*tol
                    #dφ=2pi/(3+1)
                    #φ=0:dφ:2pi-dφ
                    w=wn(ω,Γ)
                    #C=ω.+R*exp.(-1im.*φ*w) #this circle is winded opposite to Γ
                    #A+=compute_moment_matrices(L,C;l=l,N=N,output=true,random=false)
                    #moment2=compute_moment_matrices(L,C;l=l,N=N,output=true,random=false)
                    for idx=1:length(A)
                        for ll=1:Δl
                            moment=-2*pi*1.0im*w*sol.v*sol.v_adj[(idx-1)*Δl+ll]'
                            A[idx][:,ll,1]+=moment
                            A[idx][:,ll,2]+=ω*moment
                        end
                    end
                    eigvals[ω]=[sol,true]
                elseif is_new_value
                    #handle outside values
                    eigvals[ω]=[sol,false]
                end
            end

            number_of_interior_eigvals=0
            for val in values(eigvals)
                if val[2]
                    number_of_interior_eigvals+=1
                end
            end

            if neigval==number_of_interior_eigvals
                break
            else
                neigval=number_of_interior_eigvals
            end
        end
        # σ0,σ=σ,0.0
        # for idx =1:length(A)
        #     σ=max(σ,maximum(abs.(A[idx][:,:,1])))
        # end
        # σmax=max(σmax,σ)
        if loglevel>=1
            flush(stdout)
            println("maximum σmax:$σmax")
            println("final σ:$σ")
            println("relative σ/σ0:$(σ/σ0)")
            println("relative σ/σmax:$(σ/σmax)")
            println("cycles:$cycle,maxcycles:$(maxcycles<=cycle)")
        end
        #stopping criteria
        if σ/σmax<rtol_σ || σ<atol_σ
            if loglevel>=1
                println("finished after $l integrations")
            end
            break
        else
            l+=Δl
            println("you should run more outer cycles!")
        end

    end #of outer loop
    return eigvals
end
