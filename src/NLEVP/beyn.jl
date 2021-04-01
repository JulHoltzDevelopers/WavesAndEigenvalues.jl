import LinearAlgebra,FastGaussQuadrature
import ProgressMeter

"""
    Ω, P = beyn(L::LinearOperatorFamily, Γ; <keyword arguments>)

Compute all eigenvalues of `L` inside the contour `Γ` together with the associated eigenvectors. The contour is given as a list of complex numbers which are interpreted as polygon vertices in the complex plane.  The eigenvalues are stored in the list `Ω` and the eigenvectors are stored in the columns of `P`. The eigenvector `P[:,i]` corresponds to the eigenvalue `Ω[i]`, i.e., they satisfy `L(Ω[i])P[:,i]=0`.

# Arguments
- `L::LinearOperatorFamily`: Definition of the non-linear eigenvalue problem
- `Γ::Array`: List of complex points defining the contour.
- `l::Integer = 5`: estimate of the number of eigenvalues inside of `Γ`.
- `K::Integer = 1`: Augmention dimension if `Γ` is assumed to contain more than `size(L)[1]` eigenvalues.
- `N::Integer = 16`: Number of evaluation points used to perform the Gauss-Legendre integration along each edge of `Γ`.
- `tol::Float=0.0`: Threshold value to discard spurious singular values. If set to `0` (the default) no singular values are discarded.
- `pos_test::bool=true`: If set to `true` perform positions test on the computed eigenvalues, i.e., check whether the eigenvalues are enclosed by `Γ` and disregard all eigenvalues which fail the test.
- `output::bool=true`: Show progressbar if `true`.
- `random::bool=false`: Randomly initialize the V matrix if `true`.

# Returns
- `Ω::Array`: List of computed eigenvalues
- `P::Matrix`: Matrix of eigenvectors.

# Notes
The original algorithm was first presented by Beyn in [1]. The implementation closely follows the pseudocode from Buschmann et al. in [2].

# References
[1] W.-J. Beyn, An integral method for solving nonlinear eigenvalue problems, Linear Algebra and its Applications, 2012, 436(10), p.3839-3863, https://doi.org/10.1016/j.laa.2011.03.030

[2] P.E. Buschmann, G.A. Mensah, J.P. Moeck, Solution of Thermoacoustic Eigenvalue Problems with a Non-Iterative Method, J. Eng. Gas Turbines Power, Mar 2020, 142(3): 031022 (11 pages) https://doi.org/10.1115/1.4045076

See also: [`inveriter`](@ref), [`lancaster`](@ref), [`mslp`](@ref), [`rf2s`](@ref), [`traceiter`](@ref)
"""
function beyn(L::LinearOperatorFamily,Γ;l=5,K=1,N=16,tol=0.0,pos_test=true,output=true,random=false)
    #This is Beyn's algorithm as implemented in Buschmann et al. 2019
    #Beyn's original paper from 2012 also suggests a residue test which might be implemented
    #in the future
    d=size(L(0))[1]
    K=max(K,div(l,d)+Int(mod(l,d)!=0)) #ensure minimum required augmentation
    #initialize V
    if l<d #TODO consider seeding for reproducibility
        if random
            V=rand(ComplexF64,d,l)
        else
            V=zeros(ComplexF64,d,l)
            for i=1:l
                V[i,i]=1.0+0.0im
            end
        end
    else
        #TODO: there might be an easier constructor
        #V=LinearAlgebra.Diagonal(ones(ComplexF64,d))
        V=zeros(ComplexF64,d,l)
        for i=1:d
            V[i,i]=1.0+0.0im
        end
    end




    function integrand(z)
        z,w=z
        A=Array{ComplexF64}(undef,d,l,2*K)
        A[:,:,1]=L(z)\V
        A[:,:,1]*=w
        for p =1:2*K-1
            A[:,:,p+1]=z^p*A[:,:,1]
        end
        return A
    end

    A=zeros(ComplexF64,d,l,2*K) # initialize A with zeros
    A=gauss(A,integrand,Γ,N,output) #integrate over contour

    #Assemble B matrices
    B=Array{ComplexF64}(undef,d*K,l*K,2)
    rows=1:d
    cols=1:l
    for i=0:K-1,j=0:K-1
        B[rows.+d*i,cols.+l*j,1]=A[:,:,i+j+1] #indexshift +1
        B[rows.+d*i,cols.+l*j,2]=A[:,:,i+j+2]
    end

    V,Σ,W=LinearAlgebra.svd(B[:,:,1])
    #sigma check
    if output
        println("############")
        println("singular values:")
        println(Σ)
    end
    if tol>0
        mask=map(σ->σ>tol,Σ)
        V,Σ,W=V[:,mask],Σ[mask],W[:,mask]
    end
    #TODO: rank test

    #invert Σ. It's a diagonal matrix so inversion is just one line
    Σ=LinearAlgebra.Diagonal(1 ./Σ) #TODO: check zero division
    #compute eigenvalues and eigenvectors
    Ω,P = LinearAlgebra.eigen(V'*B[:,:,2]*W*Σ)
    P=V[1:d,:]*P
    #position test
    if pos_test
        mask=map(z->inpoly(z,Γ),Ω)
        Ω,P=Ω[mask],P[:,mask]
    end

    return Ω,P
end

function gauss(int,f,Γ,N,output=false)
    X,W=FastGaussQuadrature.gausslegendre(N)
    lΓ=length(Γ)
    if output
        prog = ProgressMeter.Progress(lΓ*N,desc="Beyn... ", dt=.5,
             barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=30)
    end
    for i = 1:lΓ
        if i==lΓ
            a,b=Γ[i],Γ[1]
        else
            a,b=Γ[i],Γ[i+1]
        end
        X̂=X.*(b-a)/2 .+(a+b)/2 #transform to segment
        #int+=mapreduce(f,+,zip(X̂,W)).*(b-a)/2
        iint=zero(int)
        for z in zip(X̂,W)
            iint+=f(z)
            if output
                ProgressMeter.next!(prog)
            end
        end
        int+=iint.*(b-a)/2
    end
    return int
end


#This is a Julia implementation of the point-polygon-problem as it is
#discussed on http://geomalgorithms.com/a03-_inclusion.html
#The algorithms are adapted for complex arithmetics to seamlessly
#work with beyn

function isleft(a,b,c)
    #c point to be tested
    #a, b points defining the line
    return (real(b)-real(a)) * (imag(c)-imag(a)) - (real(c)-real(a))*(imag(b)-imag(a))
end

# function inpoly(z,Γ)
#     #winding number test
#     wn=0
#     for i=1:length(Γ)
#         if i==length(Γ)
#             a,b=Γ[i],Γ[1]
#         else
#             a,b=Γ[i],Γ[i+1]
#         end
#         if imag(a) <= imag(z)
#             if imag(b)>imag(z)
#                 if isleft(a,b,z)>0
#                     wn+=1
#                 end
#             end
#         else
#             if imag(b) <= imag(z)
#                 if isleft(a,b,z)<0
#                     wn-=1
#                 end
#             end
#         end
#     end
#     return wn!=0
# end

inpoly(z,Γ)= wn(z,Γ)!=0

"""
    w=wn(z,Γ)

Compute winding number, i.e., how often `Γ` is wrapped around `z`. The sign of the winding number indocates the wrapping direction.
"""
function wn(z,Γ)
    #winding number test
    wn=0
    for i=1:length(Γ)
        if i==length(Γ)
            a,b=Γ[i],Γ[1]
        else
            a,b=Γ[i],Γ[i+1]
        end
        if imag(a) <= imag(z)
            if imag(b)>imag(z)
                if isleft(a,b,z)>0
                    wn+=1
                end
            end
        else
            if imag(b) <= imag(z)
                if isleft(a,b,z)<0
                    wn-=1
                end
            end
        end
    end
    return wn
end
## Thats Beyn in pieces...
"""
    A=compute_moment_matrices(L,Γ, <kwargs>)

Compute moment matrices of `L` integrating along `Γ`.

# Arguments
- `L::LinearOperatorFamily`
- `Γ`: List of complex points defining the contour.
- `l::Int=5`: (optional) estimate of the number of eigenvalues inside of `Γ`.
-  `K::Int=1`: (optional) Augmention dimension if `Γ` is assumed to contain more than `size(L)[1]` eigenvalues.
- `N::Int=16`: (optional)  Number of evaluation points used to perform the Gauss-Legendre integration along each edge of `Γ`.
- `output::Bool=false` (optional) toggle waitbar
- `random::Bool=false` (optional) toggle random initialization of V-matrix

# Returns
- `A::Array`: Moment matrices.

# Notes
A is a 3 dimensional array. The last index loops over the various Moments, i.e., `A[:,:,i]=∫_Γ z^i inv(L(z)) V dz`.
If `random`=false V is initialized with ones on its main diagonal. Otherwise it is random.

"""
function compute_moment_matrices(L::LinearOperatorFamily,Γ;l=5,K=1,N=16,output=false,random=false)
    d=size(L(0))[1]
    #initialize V
    if l<d #TODO consider seeding for reproducibility
        if random
            V=rand(ComplexF64,d,l)
        else
            V=zeros(ComplexF64,d,l)
            for i=1:l
                V[i,i]=1.0+0.0im
            end
        end
    else
        #TODO: there might be an easier constructor
        V=LinearAlgebra.Diagonal(ones(ComplexF64,d))
    end
    return compute_moment_matrices(L::LinearOperatorFamily,Γ,V;K=L,N=N,output=output)
end
function compute_moment_matrices(L::LinearOperatorFamily,Γ,V;K=1,N=16,output=false)
    d,l=size(V)

    function integrand(z)
        z,w=z
        A=Array{ComplexF64}(undef,d,l,2*K)
        A[:,:,1]=L(z)\V
        A[:,:,1]*=w
        for p =1:2*K-1
            A[:,:,p+1]=z^p*A[:,:,1]
        end
        return A
    end
    A=zeros(ComplexF64,d,l,2*K) # initialize A with zeros
    A=gauss(A,integrand,Γ,N,output) #integrate over contour

    return A
end
"""
    Ω,P = moments2eigs(A; tol_σ=0.0, return_σ=false)

Compute eigenvalues `Ω` and eigenvectors `P from moment matrices `A`.

# Arguments
- `A::Array`: moment matrices
- `tol::Float=0.0`: tolerance for truncating singular values.
- `return_σ::Bool=false`: if true also return the singular values as third return value.

# Returns
- `Ω::Array`: List of computed eigenvalues
- `P::Matrix`: Matrix of eigenvectors.
- `Σ::Array`: List of singular values (only if `return_σ==true`)

# Notes
The the i-th coloumn in P corresponds to the i-th eigenvalue in `Ω`, i.e., `Ω[i]` and `P[:,i]` are an eigenpair.

See also: [`compute_moment_matrices`](@ref)
"""
function moments2eigs(A;tol_σ=0.0,return_σ=false)
    #Assemble B matrices
    d=size(A[1],1)
    Δl=size(A[1],2)
    l=length(A)*Δl
    K=size(A[1],3)÷2
    B=Array{ComplexF64}(undef,d*K,l*K,2)
    rows=1:d
    cols=1:Δl
    for i=0:K-1,j=0:K-1
        for ll=1:length(A)
            B[rows.+d*i,cols.+(ll-1).*Δl.+l.*j,1].=A[ll][:,:,i+j+1] #indexshift +1
            B[rows.+d*i,cols.+(ll-1).*Δl.+l.*j,2].=A[ll][:,:,i+j+2]
        end
    end
    V,Σ,W=LinearAlgebra.svd(B[:,:,1])
    #sigma check
    if tol_σ>0
        mask=map(σ->σ>tol,Σ)
        V,Σ,W=V[:,mask],Σ[mask],W[:,mask]
    end
    #println("######")
    #println(Σ)
    #invert Σ. It's a diagonal matrix so inversion is just one line
    invΣ=LinearAlgebra.Diagonal(1 ./Σ) #TODO: check zero division
    #compute eigenvalues and eigenvectors
    Ω,P = LinearAlgebra.eigen(V'*B[:,:,2]*W*invΣ)
    P=V[1:d,:]*P

    if return_σ
        return Ω,P,Σ
    else
        return Ω,P
    end
end


"""
    Ω,P=pos_test(Ω,P,Γ)

Filter eigenpairs `Ω` and `P`for those whose eigenvalues lie inside `Γ`.

See also: [`moments2eigs`](@ref)
"""
function pos_test(Ω,P,Γ)
    mask=map(z->inpoly(z,Γ),Ω)
    Ω,P=Ω[mask],P[:,mask]
    return Ω,P
end

## count poles and zeros
"""
    n=count_poles_and_zeros(L,Γ;N=16,output=false)

Count number "n" of poles and zeros of the determinant of `L` inside the contour `Γ`.
The optional parameters `N` and `output`. respectively  determine the order of the
Gauss-Legendre integration and toggle whether output is desplayed.

# Notes
Pole orders count negative while zero order count positive.
The evaluation is based on applying the residue theorem to the determinant.
It utilizes Jacobi's formula to compute the necessary derivatives from a
trace operation. and is therefore only suitable for small problems.


"""
function count_poles_and_zeros(L,Γ;N=16,output=false)
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
## Projection method stuff
"""
    V=initialize_V(d::Int,l::Int,random::Bool=false)

initialioze matrix with `d`rows and `l` colums either as diagonal matrix
featuring ones on the main diagonal and anywhere else or with random entries
(`random`=true).

See also: [`generate_subspace`](@ref)
"""
function initialize_V(d::Int,l::Int;random::Bool=false)
    if random
        V=rand(ComplexF64,d,l)
        for i=1:l
            V[:,i].\=LinearAlgebra.norm(V[:,i])
        end
    else
        V=zeros(ComplexF64,d,l)
        for i=1:l
            V[i,i]=1.0+0.0im
        end
    end
    return V
end
raw"""
    Q,resnorm=generate_subspace(L,Y,tol,Z,N;output=true))

Compute orthonormal `Q` basis for a subspace that can be used to reduce the
dimension of `L`. The created subspace gurantees the residual of `L(z)\V` to be less than
`tol` for each `z∈Z`.

# Arguments
- `L::LinearOperatorFamily`: operator family for which the subspace is to be computed
- `Y::matrix`: matrix for residual test
- `tol::Float64`: tolerance for residual test
- `Z::List`: List of sample points
- `N::Int`: (optional) if set the list Z is interpreted as edges of a closed polygonal contour and N Gauss-Legendre sample points are generated for each edge.
- `output::Bool=false`: (optional) toggle progressbar
- `include_Y::Bool=true`: (optional) include Y in subspace

# Returns

- `Q::Matrix`: unitary matrix whose colums form the basis of the subspace
- `resnorm::List`: list of residuals at the sample points computed during subspace generation. (these values are an upper bound)

# Notes

The algorithm is based on the idea in [1]. However, it uses incremental QR
decompositions to built the subspace and is not computing Beyn's integral.
However, the obtained matrix can be used to project `L` on teh subspace and use
Beyn's algorithm or any other eigenvalue solver on the projected problem.

# References

[1] A Study on Matrix Derivatives for the Sensitivity Analysis of Electromagnetic
Eigenvalue Problems, P. Jorkowski and R. Schuhmann, 2020, IEEE Trans. Magn., 56,
[doi:10.1109/TMAG.2019.2950513](https://doi.org/10.1109/TMAG.2019.2950513)

See also: [`beyn`](@ref), [`initialize_V`](@ref), [`project`](@ref)
"""
function generate_subspace(L,Y,tol,Z;output::Bool=true,tol_err::Float64=Inf,include_Y=true)
    #TODO: sanity checks for sizes
    d,k=size(Y)
    dim=k
    N=length(Z) #IDEA: consider random permutation
    ## initialize subspace from first point
    A=[] #compressed qr storage for incremental qr
    τ=[] #compressed qr storage for incremental qr

    for kk = 1: k
        if include_Y
            qrfactUnblockedIncremental!(τ,A,Y[:,kk:kk])
        else
            qrfactUnblockedIncremental!(τ,A,L(Z[1])\Y[:,kk:kk])
        end
    end
    #QR=LinearAlgebra.qr(L(Z[1])\Y)
    #Q=QR.Q[:,1:dim]# TODO: optimize these QR calls
    resnorm=zeros(N*k)
    idx=0
    Q=get_Q(τ,A) #TODO: incrementally create only the last column, consider preallocation for speed
    QY=Q'*Y #project rhs


    ## Preallacations
    #Y=Array{ComplexF64}(undef,d,1)
    #Y_exact=Array{ComplexF64}(undef,d,1)

    if output
        prog = ProgressMeter.Progress(k*N,desc="Subspace.. ", dt=1,
         barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),)
    end
    for z in Z
        if dim==d
            break
        end
        idx+=1
        Lz=L(z)
        #Lnorm=LinearAlgebra.norm(Lz)
        Lnorm=1
        QLQ=Q'*Lz*Q #projection to subspace
        for kk=1:k
            x=QLQ\QY[:,kk] #solve projected problem
            X=Q*x #backprojection
            #residual test
            #weights=ones(ComplexF64,d)
            #weights[20:end].=0
            #weights=get_weights(Lz)
            weights=1
            res=LinearAlgebra.norm((weights).*(Lz*X-Y[:,kk]))
            res/=Lnorm
            #println("$(kk+(idx-1)*k): $res vs $tol")
            #TODO: check that supdspace does not get linearly dependent
            # in order to enable 0 tolerance
            if res>tol
                # if true# idx<=3
                #     println("vorher:")
                #     println(idx," ",kk+(idx-1)*k," ",res," ",resnorm[5])
                #     println("nacher:")
                #     flush(stdout)
                # end
                dim+=1
                X_exact=Lz\Y[:,kk:kk]
                #err=LinearAlgebra.norm(weights.*(X_exact.-X))
                #println("error:$err")
                X=X_exact
                res=LinearAlgebra.norm((weights).*(Lz*X-Y[:,kk]))
                #if err>tol_err
                #    #tol=max(tol,res/10) #TODO: the minimum operation is not necessary
                #    tol=res*10
                #end

                #res=err
                res/=Lnorm
                #IMPORTANT qrfactUnblockedIncremental! overwrites X!
                #Therfore command must come after residual calculation
                qrfactUnblockedIncremental!(τ,A,X[:,1:1])
                #Q=create_Q(τ,A)
                Q_old,Q=Q,Array{ComplexF64}(undef,d,dim)
                Q[:,1:end-1]=Q_old[:,:]
                #initialize last column with e_dim_vector
                #this is unused but allacated space. Therefore its safe to be used
                # e_dim will be overwritten later to with the true
                # column, but it is used to initialize finding the
                #true column.
                Q[:,end]=zeros(ComplexF64,d)
                Q[dim,end]=1
                Q[:,end]=mult_QV(τ,A,Q[:,end]) #build last column
                ##reinitialize cache
                QLQ=Q'*Lz*Q #projection to subspace #TODO: do this incrementally using last column only
                QY=Q'*Y# project rhs
                #ΔQLQ=Q*L*ΔQ #only build last colum



                #QV=[QV; ΔQ'*V]
                # if true#idx<=3
                #     println(idx," ",kk+(idx-1)*k," ",res," ",resnorm[5])
                #     flush(stdout)
                # end
            end
            resnorm[kk+(idx-1)*k]=res
            if output
                ProgressMeter.next!(prog)
            end
        end
    end
    #Q=create_Q(τ,A)
    if output
        println("Finished subspace generation!")
        if tol_err<Inf
            println("adapted residual tolerance is $tol .")
        end
    end
    return Q,resnorm
end

function generate_subspace(L,Y,tol,Γ,N::Int;output::Bool=true,tol_err::Float64=Inf,include_Y=true)
    ## generate list of Gauss-Legendre points
    X,W=FastGaussQuadrature.gausslegendre(N)
    lΓ=length(Γ)
    Z=zeros(ComplexF64, lΓ*N) #TODO: undef intialization
    #TODO: remove doubled edge points!
    for i = 1:lΓ
        if i==lΓ
            a,b=Γ[i],Γ[1]
        else
            a,b=Γ[i],Γ[i+1]
        end
        Z[1+(i-1)*N:i*N]=X.*(b-a)/2 .+(a+b)/2 #transform to segment
    end
    return generate_subspace(L,Y,tol,Z,output=output,tol_err=tol_err,include_Y=include_Y)
end


"""
    P::LinearOperatorFamily=project(L::LinearOperatorFamily,Q)

Project `L` on the subspace spanned by the unitary matrix `Q`, i.e.
`P(z)=Q'*L(z)*Q`.

See also: [`generate_subspace`](@ref)
"""
function project(L::LinearOperatorFamily,Q)
    P=LinearOperatorFamily()
    P.params=deepcopy(L.params)
    P.eigval=L.eigval
    P.auxval=L.auxval
    P.mode=deepcopy(L.mode)
    P.active=deepcopy(L.active)

    #term=Term(A2,(pow2,),((:λ,),),"A2")
    for term in L.terms
        M=Q'*term.coeff*Q
        push!(P,Term(M,deepcopy(term.func),deepcopy(term.params),
            deepcopy(term.symbol),deepcopy(term.operator)))
    end
    return P
end

## incremental QR
# TODO: merge request for julia base
#===
hack of line 185 in
https://github.com/JuliaLang/julia/blob/69fcb5745bda8a5588c089f7b65831787cffc366/stdlib/LinearAlgebra/src/qr.jl#L301-L378
for obtaining incremental qr
===#
using LinearAlgebra: reflector!, reflectorApply!
##
function qrfactUnblockedIncremental!(τ,A,V::AbstractMatrix{T}) where {T} #AbstarctArray?
    #require_one_based_indexing(A)
    #require_one_based_indexing(T)
    V=deepcopy(V)
    m=length(A)
    n=length(V)
    #apply all previous reflectors
    for k=1:m
        x = view(A[k], k:n,1)#TODO: consider end
        reflectorApply!(x,τ[k],view(V,k:n,1:1))
    end
    #create new reflector
    x = view(V, m+1:n,1)
    τk= reflector!(x)

    #append latest updates

    #if  !isapprox(V[m+1],0.0+0.0im,atol=1E-8)
        push!(τ,τk)
        push!(A,V)
    #end
    return
end
##
function mult_VQ(τ,A,V)
    m=length(A)
    n=length(A[1])
    #@inbounds begin
        for k=1:m
            τk=τ[k]
            vk=zeros(ComplexF64,n)
            vk[k]=1
            vk[k+1:end]=A[k][k+1:end]
            V=(V-(V*vk)*(vk'*τk))
        end
    #end
    return V
end
function mult_QV(τ,A,V;trans=true)
    m=length(A)
    n=length(A[1])
    M=1:m
    if trans
        M=reverse(M)
    end
    #@inbounds begin
        for k=M #its transposed so householder transforms are executed in reverse order
            τk=τ[k]
            vk=zeros(ComplexF64,n)
            vk[k]=1
            vk[k+1:end]=A[k][k+1:end]
            V=(V-(τk*vk)*(vk'*V))
        end
    #end
    return V
end

function get_Q(τ,A)
    n=size(A[1],1)
    dim=length(τ)
    mult_QV(τ,A,initialize_V(n,dim))
end

function get_R(A)
    n=length(A)
    m=length(A[1])
    R=zeros(ComplexF64,m,n)
    for i=1:n
        R[1:i,i]=A[i][1:i]
    end
    return R
end


function get_weights(M)
    m,n=size(M)
    weights=Array{ComplexF64}(undef,m)
    for i=1:m
        weights[i]=M[i,i]#sum(M[i,:])
    end
    return 1. /weights
end


## Givens rotation
function givens_rot!(M,i,j)
    a=M[i]
    b=M[j]
    if b!=0
        #r=hypot(a,b)
        #hack of chypot
        r = a/b
        r = b*sqrt(1.0+0.0im+r*r)
        c=a/r
        s=-b/r
    else
        c = 1.0+0.0im
        s = 0.0+0.0im
        r = a
    end
    ϕ=log(c+1.0im*s)/1.0im
    #apply givens rotation
    M[i]=r
    M[j]=ϕ
    return nothing
end
function incremental_QR!(A,v)
    push!(A,v)
    n=length(A)
    m=length(A[1])
    #apply previous givens rotations
    for j=1:n-1
        for i=j+1:m
            ϕ=A[j][i]
            c=cos(ϕ)
            s=sin(ϕ)
            A[n][j],A[n][i]=c*A[n][j]-s*A[n][i], s*A[n][j]+c*A[n][i]
        end
    end

    #do new givens rotations
    for j=n+1:m
        givens_rot!(A[n],n,j)
    end
    return nothing
end
function get_Q(A)
    m=length(A[1])
    n=length(A)
    Q=zeros(ComplexF64,m,n)
    for i=1:n
        Q[i,i]=1
    end
    for i=1:m
        for j=i+1:n
            ϕ=A[i][j]
            c=cos(ϕ)
            s=sin(ϕ)
            Q[:,i],Q[:,j]=c*Q[:,i]-s*Q[:,j], s*Q[:,i]+c*Q[:,j]
        end
    end
    return Q
end


#TODO: hack the imporved julia version
function chypot(a,b)
    if a == 0
        h=0
    else
        r = b/a
        h = a*sqrt(1+r*r)
    end
    return h
end
