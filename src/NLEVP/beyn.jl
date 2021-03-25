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
