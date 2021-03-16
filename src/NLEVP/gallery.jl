module Gallery
import LinearAlgebra
using ..NLEVP

"""
    D,x=cheb(N)

Compute differentiation matrix `D` on Chebyshev grid `x`. See Lloyd N.
Trefethen. Spectral methods in MATLAB. Society for Industrial and Applied
Mathematics, Philadelphia, PA, USA, 2000.
"""
function cheb(N)
    # Compute D = differentiation matrix, x = Chebyshev grid.
    # See Lloyd N. Trefethen. Spectral Methods in MATLAB. Society for Industrial
    # and Applied Mathematics, Philadelphia, PA, USA, 2000.
    if N == 0
        D=0
        x=1
        return D,x
    end

    x = cos.(pi/N*(0:N))
    c = [2.; ones(N-1); 2. ].*((-1.).^(0:N))
    X = repeat(x,1,N+1)
    dX= X-X'
    I=Array{Float64}(LinearAlgebra.I,N+1,N+1)
    D =(c*(1. ./c'))./(dX+I) # off-diagonal entries
    for i = 1:size(D)[1]
        D[i,i]-=sum(D[i,:])
    end
    #D=D -np.diag(np.array(np.sum(D,1))[:,0]) #diagonal entries
    return D,x
end

"""
    L,y=orr_sommerfeld(N=64, Re=5772, w=0.26943)

Return a linear operator family `L` representing the Orr-Sommerfeld equation.
A Chebyshev collocation method is used to discretize the problem. The grid
points are returned as `y`. The number of grid points is `N`.

#Note
The OrrSommerfeld equation reads:
[(d²/dy²-λ²)²-Re((λU-ω)(d²/dy²-λ²)-λU'')] v = 0

where y is the spatial coordinate, λ the wavenumber (set as the eigenvalue by
default), Re the Reynolds number, ω the oscillation frequency, U the mean flow
velocity profile and v the velocity fluctuation amplitude in the y-direction
(the mode shape).

#Remarks on implementation
This is a Julia port of (our python port of) the Orr-Sommerfeld example from the
Matlab toolbox "NLEVP: A Collection of Nonlinear Eigenvalue Problems"
by  T. Betcke, N. J. Higham, V. Mehrmann, C. Schröder, and F. Tisseur.

The original toolbox is available at :
http://www.maths.manchester.ac.uk/our-research/research-groups/numerical-analysis-and-scientific-computing/numerical-analysis/software/nlevp/

See  [1] or [2] for further refrence.

#References
[1] T. Betcke, N. J. Higham, V. Mehrmann, C. Schröder, and F. Tisseur, NLEVP:
A Collection of Nonlinear Eigenvalue Problems, MIMS EPrint 2011.116, December
2011.

[2] T. Betcke, N. J. Higham, V. Mehrmann, C. Schröder, and F. Tisseur, NLEVP:
A Collection of Nonlinear Eigenvalue Problems. Users' Guide, MIMS EPrint
2010.117, December 2011.
"""
function orr_sommerfeld(N=64, Re=5772, ω =0.26943)
    # Define the Orr-Sommerfeld operator for spatial stability analysis.
    # N: number of Cheb points, R: Reynolds number, w: angular frequency
    N=N+1
    # 2nd- and 4th-order differentiation matrices:
    D,y=cheb(N); D2 =D^2;  D2 = D2[2:N,2:N]
    S= LinearAlgebra.diagm(0=>[0; 1.0./(1. .-y[2:N].^2.); 0])
    D4=(LinearAlgebra.diagm(0=>1. .-y.^2.)*D^4-8*LinearAlgebra.diagm(0=>y)*D^3-12*D^2)*S
    D4 = D4[2:N,2:N]

    I=Array{ComplexF64}(LinearAlgebra.I,N-1,N-1)
    # type conversion
    D2=D2.+0.0im
    D4=D4.+0.0im
    U=LinearAlgebra.diagm(0=>-y[2:N].^2.0.+1)
    # Build Orr-Sommerfeld operator
    L= LinearOperatorFamily(["λ","ω","Re","a"],complex([1.,ω,Re,Inf]))
    push!(L,Term(I,(pow_a(4),),((:λ,),),"λ^4","I"))
    push!(L,Term(1.0im*U,(pow_a(3),pow1,),((:λ,),(:Re,),),"iλ^3Re","i*U"))
    push!(L,Term(-2*D2,(pow2,),((:λ,),),"λ^2","-2D2"))
    push!(L,Term(-1.0im*I,(pow2,pow1,pow1,),((:λ,),(:ω,),(:Re,),),"λ^2*ω*Re","-i*I"))
    push!(L,Term(-1.0im*(U*D2+2.0*I),(pow1,pow1,),((:λ,),(:Re,),),"λ*Re","(U*D2+2*I)"))
    push!(L,Term(1.0im*D2,(pow1,pow1),((:ω,),(:Re,),),"ω*Re","i*D2"))
    push!(L,Term(D4,(),(),"","D4"))
    push!(L,Term(-I,(pow1,),((:a,),),"-a","__aux__"))
    return L,y
end

"""
    L,x,y=biharmonic(N=12;scaleX=2,scaleY=1+sqrt(5))

Discretize the biharmonic equation.

# Arguments
-`N::Int`: Number of collocation points
-`scaleX::Float=2` length of the x-axis
-`scaleY::Float=1+sqrt(5)` length of the y-axis

# Returns
-`L::LinearOperatorFamily`: discretization of the biharmonic operator
- `x:Array`: x-axis
- `y:Array`: y-axis

# Notes
The biharmonic equation models the oscillations of a membrane. It is an an eigenvalue problem that reads:
(∇⁴+εcos(2πx)cos(πy))**v**=λ**v** in Ω
with boundary conditions
**v**=∇²**v**=0


The term εcos(2πx)cos(πy) is included to model some inhomogeneous material properties.
The equation is discretized using Chebyshev collocation.

See also: [`cheb`](@ref)
"""

function biharmonic(N=12;scaleX=2,scaleY=1+sqrt(5))
#N,scaleX,scaleY=12,2,1+sqrt(5)
    N=N+1
    D, xx = cheb(N)
    x = xx/scaleX
    y = xx/scaleY

    #The Chebychev matrices are space dependent scale the accordingly
    Dx = D*scaleX
    Dy = D*scaleY

    D2x=Dx*Dx
    D2y=Dy*Dy
    Dx = Dx[2:N,2:N]
    Dy = Dy[2:N,2:N]
    D2x= D2x[2:N,2:N]
    D2y= D2y[2:N,2:N]
    #Apply BC
    I = Array{ComplexF64}(LinearAlgebra.I,N-1,N-1)
    L = kron(I,D2x) +kron(D2y,I)
    X=kron(ones(N-1),x[2:N])
    Y=kron(y[2:N],ones(N-1))
    P=LinearAlgebra.diagm(0=>cos.(π*2*X).*cos.(π*Y))
    D4= L*L
    I = Array{ComplexF64}(LinearAlgebra.I,(N-1)^2,(N-1)^2)
    L= LinearOperatorFamily(["λ","ε","a"],complex([0,0.,Inf]))
    push!(L,Term(D4,(),(),"","D4"))
    push!(L,Term(P,(pow1,),((:ε,),),"ε","P"))
    push!(L,Term(-I,(pow1,),((:λ,),),"-λ","I"))
    push!(L,Term(-I,(pow1,),((:a,),),"-a","__aux__"))
    return L,x,y
end


## 1DRijkeModel
using SparseArrays
"""
    L,grid=rijke_tube(resolution=128; l=1, c_max=2,mid=0)

Discretize a 1-dimensional Rijke tube. The model is based on the thermoacoustic equations.
This eigenvalue problem reads:
∇c²(x)∇p+ω²p-n exp(-iωτ) ∇p(x_ref)=0 for x in ]0,l[
with boundariy conditions
∇p(0)=p(l)=0
"""
function rijke_tube(resolution=127; l=1, c_max=2,mid=0)
  n=1.0
  tau=2
  #l=1
  c_min=1
  #c_max=2
  outlet=resolution
  outlet_c=c_max
  grid=range(0,stop=l,length=resolution)
  e2p=[(i, i+1) for i =1:resolution-1]
  number_of_elements=resolution-1
  if mid==0
    mid=div(resolution,2)+1# this is the element containing the flame
    #it is not the center if resolotion is odd but then the refernce is in the center
  end
  ref=mid-1 #the reference element is the one in the middle
  flame=(mid)
  #compute element volume
  e2v=diff(grid)
  V=e2v[mid]
  e2c=[i < mid ? c_min : c_max for i = 1:resolution]

  #assemble mass matrix
  m_unit=[2 1;1 2]*1/6
  #preallocation
  ii=Array{Int64}(undef,(resolution-1,2^2))
  jj=Array{Int64}(undef,(resolution-1,2^2))
  mm=Array{ComplexF64}(undef,(resolution-1,2^2))
  for (idx,el) in enumerate(e2p)
    mm[idx,:] = m_unit[:]*e2v[idx]
    ii[idx,:] = [el[1] el[2] el[1] el[2]]
    jj[idx,:] = [el[1] el[1] el[2] el[2]]
  end
  # finally assemble global (sparse) mass matrix
  M =sparse(ii[:],jj[:],mm[:])

  # assemble stiffness matrix
  k_unit = -[1. -1.; -1. 1.]
  #preallocation
  ii=Array{Int64}(undef,(resolution-1,2^2))
  jj=Array{Int64}(undef,(resolution-1,2^2))
  kk=Array{ComplexF64}(undef,(resolution-1,2^2))
  for (idx,el) in enumerate(e2p)
    kk[idx,:]=k_unit[:]/e2v[idx]*e2c[idx]^2
    ii[idx,:] = [el[1] el[2] el[1] el[2]]
    jj[idx,:] = [el[1] el[1] el[2] el[2]]
  end
  # finally assemble global (sparse) stiffness matrix
  K = sparse(ii[:],jj[:],kk[:])


  #assemble boundary mass matrix
  bb=[-outlet_c*1im]
  ii=[outlet]
  jj=[outlet]
  B = sparse(ii,jj,bb)
  #assemble flame matrix
  #preallocation
  ii=Array{Int64}(undef,(length(flame),2*2))
  jj=Array{Int64}(undef,(length(flame),2*2))
  qq=Array{ComplexF64}(undef,(length(flame),2*2))

  grad_p_ref=[-1 1]
  grad_p_ref=grad_p_ref/e2v[ref]*1
  #println("##!!!!!Here!!!!!##'")
  #println("flame:",flame)
  #println("ref:",ref)
  for (idx,el) in enumerate(flame)
    jj[idx,:]=[e2p[ref][1] e2p[ref][2] e2p[ref][1] e2p[ref][2]]
    ii[idx,:]=[e2p[el][1] e2p[el][1] e2p[el][2] e2p[el][2]]
    run=1
    for i=1:2
      for j =1:2
        qq[idx,run]= grad_p_ref[1,j]*e2v[el]/2
        run+=1
      end
    end
  end
  #finally assemble flame matrix
  Q= sparse(ii[:],jj[:],-qq[:],resolution,resolution)
  Q=Q/V
  L=LinearOperatorFamily(["ω","n","τ","Y","λ"],complex([0.,n,tau,1E15,Inf]))
  push!(L,Term(M,(pow2,),((:ω,),),"ω^2","M"))
  push!(L,Term(K,(),(),"","K"))
  push!(L,Term(B,(pow1,pow1,),((:ω,),(:Y,),),"ω*Y","C"))
  push!(L,Term(Q,(pow1,exp_delay,),((:n,),(:ω,:τ,)),"n*exp(-i ω τ)","Q"))
  push!(L,Term(-M,(pow1,),((:λ,),),"-λ","__aux__"))
  #push!(L,Term(-SparseArrays.I,(pow1,),((:λ,),),"-λ","__aux__"))
 return L,grid
end

end
