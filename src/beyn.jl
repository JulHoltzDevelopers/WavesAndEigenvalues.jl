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
- `output::bool=false`: Show progressbar if `true`.
# Returns
- `Ω::Array`: List of computed eigenvalues
- `P::Matrix`: Matrix of eigenvectors.

# Notes
The original algorithm was first presented by Beyn in [1]. The implementation closely follows the pseudocode from Buschmann et al. in [2].

# References
[1] W.-J. Beyn, An integral method for solving nonlinear eigenvalue problems, Linear Algebra and its Applications, 2012, 436(10), p.3839-3863, https://doi.org/10.1016/j.laa.2011.03.030

[2] P.E. Buschmann, G.A. Mensah, J.P. Moeck, Solution of Thermoacoustic Eigenvalue Problems with a Non-Iterative Method, J. Eng. Gas Turbines Power, Mar 2020, 142(3): 031022 (11 pages) https://doi.org/10.1115/1.4045076

See also: [`householder`](@ref), [`juniper`](@ref), [`nicoud`](@ref), [`picard`](@ref)
"""
function beyn(L::LinearOperatorFamily,Γ;l=5,K=1,N=16,tol=0.0,pos_test=true,output=false)
    #This is Beyn's algorithm as implemented in Buschmann et al. 2019
    #Beyn's original paper from 2012 also suggests a residue test which might be implemented
    #in the future
    d=size(L(0))[1]
    #initialize V
    if l<d #TODO consider seeding for reproducibility
        V=rand(ComplexF64,d,l)
    else
        #TODO: there might be an easier constructor
        V=LinearAlgebra.Diagonal(ones(ComplexF64,d))
    end

    function integrand(z)
        z,w=z
        #A=[Array{ComplexF64}(undef,dims) for p=0:2*K-1]
        A=zeros(ComplexF64,d,l,2*K,) #TODO no zero initialization
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
    if tol>0
        mask=map(σ->σ>tol,Σ)
        V,Σ,W=V[:,mask],Σ[mask],W[:,mask]
    end
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
        p = ProgressMeter.Progress(lΓ*N,desc="Beyn... ", dt=.5,
             barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=10)
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
                ProgressMeter.next!(p)
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

function inpoly(z,Γ)
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
    return wn!=0
end
