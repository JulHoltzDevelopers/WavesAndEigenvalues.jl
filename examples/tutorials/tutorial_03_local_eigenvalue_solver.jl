# #Tutorial 03
#
# ##Introduction
# This tutorial will teach you the details of the local eigenvalue solver.
# The solver is a generalization of Lancasters *Generalised Rayleigh Quotient
# iteration* [1], in the sense that it enables not only Newton's method for root
# finding, but also up to fith-order Householder iterations. The implementation
# closely follows the consideration in [2].
#
#[1] P. Lancaster, A Generalised Rayleigh Quotient Iteration for Lambda-Matrices,Arch. Rational Mech Anal., 1961, 8, p.
#309-322, https://doi.org/10.1007/BF00277446
#
#[2] G.A. Mensah, Efficient Computation of Thermoacoustic Modes, Ph.D. Thesis, TU Berlin, 2019
## #jl
# ##Model set up.
# The model is the same Rijke tube configuration as in Tutorial 01:
using WavesAndEigenvalues.Helmholtz
mesh=Mesh("Rijke_mm.msh",scale=0.001) #load mesh
dscrp=Dict() #initialize model discreptor
dscrp["Interior"]=(:interior, ()) #define resonant cavity
dscrp["Outlet"]=(:admittance, (:Y,1E15)) #specify outlet BC
γ=1.4 #ratio of specific heats
ρ=1.225 #density at the reference location upstream to the flame in kg/m^3
Tu=300.0    #K unburnt gas temperature
Tb=1200.0    #K burnt gas temperature
P0=101325.0 # ambient pressure in Pa
A=pi*0.025^2 # cross sectional area of the tube
Q02U0=P0*(Tb/Tu-1)*A*γ/(γ-1) #the ratio of mean heat release to mean velocity Q02U0
x_ref=[0.0; 0.0; -0.00101] #reference point
n_ref=[0.0; 0.0; 1.00] #directional unit vector of reference velocity
n=0.01 #interaction index
τ=0.001 #time delay
dscrp["Flame"]=(:flame,(γ,ρ,Q02U0,x_ref,n_ref,:n,:τ,n,τ)) #flame dynamics
R=287.05 # J/(kg*K) specific gas constant (air)
speedofsound(x,y,z) = z<0. ? sqrt(γ*R*Tu) : sqrt(γ*R*Tb)
c=generate_field(mesh,speedofsound)
L=discretize(mesh,dscrp,c)

# ## The local solver
#
# The key idea behind the local solver is simple: The eigenvalue problem
# `L(ω)p==0` is reinterpreted as `L(ω)p == λ*Y*p`. This means, the nonlinear
# eigenvalue `ω` becomes a parameter in a linear eigenvalue problem with
# (auxiliary) eigenvalue `λ`. If an `ω` is found such that `λ==0`, this very
# `ω` solves the original non-linear eigenvalue problem. Because the auxiliary
# eigenvalue problem is linear the implicit relation `λ=λ(ω)` can be examined
# with established perturbation theory of linear eigenvalue problems. This
# theory can be used to compute the derivative `dλ/dω` and iteratively find the
# root of `λ=λ(ω)`. Each of these iterations requires solving a linear eigenvalue
# problem and higher order derivatives improve the  iteration. The options of
# the local-solver `householder` allow the user to control the solution process.
#
# ## Mandatory inputs
#
# Mandatory input arguments are only the linear operator family `L` and an
# initial guess `z` for the eigenvalue `ω`` iteration.
z=300.0*2*pi
sol,nn,flag = householder(L, z)

# ## Maximum iterations
#
# As per default five iteration are performed. This iteration number can be
# customized using the keyword `maxiter` For instance the following command
# runs 30 iterations
sol,nn,flag = householder(L, z, maxiter=30)

# Terminating a convergend solutions
#
# Usually, the procedure converges after a few iteration steps and further
# iterations do not significantly improve the solution. For this reason the
# keyword `tol` allows to set a threshold, such that two consecutive iterates
# for the eigenvalue fullfill `abs(ω_0-ω_1)`, the iteration is terminated.
# This keyword defaults to `tol=0.0` and, thus, `maxiter`iterations will be
# performed. However, it is highly recommended to set the parameter to some
# positive value whenever possible.
sol,nn,flag = householder(L, z, maxiter=30, tol=1E-10)
# Note that the toloerance is also a bound for the error associated with the
# computed eigenvalue. Tip: The 64-bit floating numbers have an accuracy of
# `eps≈1E-16`. Hence, the threshold `tol` should not be chosen less than 16
# orders of magnitude smaller than the expected eigenvalue. Indeed, `tol=1E-10`
# is already close to the machine-precision we can expect for the Rijke-tube
# model.
#
# ##Determining convergence
#
# Technically, the iteration may stop because of slow progress in the itreation
# rather than actual convergence. A simple indicator for actual convergnce is
# the auxiliary eigenvalue `λ`. The closer it is to `0` the better is the
# quality of the computed solution. To help the identification of falsely
# terminated iterations you can specify another tolerance `lam_tol`. If
# `abs(λ)>lam_tol` the computed eigenvalue is deemed to be spurious. This
# feature only makes sense when a termination threshold has been specified.
# As in the following example:
sol,nn,flag = householder(L, z, maxiter=30, tol=1E-10, lam_tol=1E-8)
# ##Faster convergence
# In order to improve the convergence rate, higher-order derivatives may be
# computed. You can use up to fifth order perturbation theory to improve the
# iteration via the `order` keyword. For instance, third-order theory is used
# in this example
sol,nn,flag = householder(L, z, maxiter=30, tol=1E-10, lam_tol=1E-8, order=3)
#
# ##
#
# Under the hood the routine calls ARPACK to utilize Arnoldi's method to solve
# the linear (auxiliary) eigenvalue problem. This method can compute more than
# one eigenvalue close to some initial guess. Per default, only one is sought
# but via the `n_eig_val` keyword the number can be increased. This results
# in an increased computation time but provides more canditates for the next
# iteration, potentialy improving the convergence.
sol,nn,flag = householder(L, z, maxiter=30, tol=1E-10, lam_tol=1E-8, n_eig_val=3)
#
#TODO: explain solution-object and error flags, degeneracy, relaxation, v0
