```@meta
EditURL = "<unknown>/tutorial_03_local_eigenvalue_solver.jl"
```

# Tutorial 03  A Local Eigenvalue Solver

This tutorial will teach you the details of the local eigenvalue solver.
The solver is a generalization of the method of successive linear problems
[1], in the sense that it enables not only Newton's method for root
finding, but also high order Householder iterations. The implementation
closely follows the considerations in [2].

## Model set-up.
The model is the same Rijke tube configuration as in Tutorial 01:

```julia
using WavesAndEigenvalues.Helmholtz
mesh=Mesh("Rijke_mm.msh",scale=0.001) #load mesh
dscrp=Dict() #initialize model descriptor
dscrp["Interior"]=(:interior, ()) #define resonant cavity
dscrp["Outlet"]=(:admittance, (:Y,1E15)) #specify outlet BC
γ=1.4 #ratio of specific heats
ρ=1.225 #density at the reference location upstream to the flame in kg/m^3
Tu=300.0    #K unburnt gas temperature
Tb=1200.0    #K burnt gas temperature
P0=101325.0 #ambient pressure in Pa
A=pi*0.025^2 #cross sectional area of the tube
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
```

```
1006×1006-dimensional operator family: 

ω^2*M+K+n*exp(-iωτ)*Q+ω*Y*C

Parameters
----------
n	0.01 + 0.0im
λ	Inf + 0.0im
ω	0.0 + 0.0im
τ	0.001 + 0.0im
Y	1.0e15 + 0.0im

```

## The local solver

The key idea behind the local solver is simple: The eigenvalue problem
`L(ω)p==0` is reinterpreted as `L(ω)p == λ*Y*p`. This means, the nonlinear
eigenvalue `ω` becomes a parameter in a linear eigenvalue problem with
(auxiliary) eigenvalue `λ`. If an `ω` is found such that `λ==0`, this very
`ω` solves the original non-linear eigenvalue problem. Because the auxiliary
eigenvalue problem is linear the implicit relation `λ=λ(ω)` can be examined
with established perturbation theory of linear eigenvalue problems. This
theory can be used to compute the derivative `dλ/dω` and iteratively find the
root of `λ=λ(ω)`. Each of these iterations requires solving a linear eigenvalue
problem and higher order derivatives improve the  iteration. The options of
the local-solver `mslp` allows the user to control the solution process.

## Mandatory inputs

Mandatory input arguments are only the linear operator family `L` and an
initial guess `z` for the eigenvalue `ω`` iteration.

```julia
z=300.0*2*pi
sol,nn,flag = mslp(L, z)
```

```
(####Solution####
eigval:
ω = 1710.69777723934 + 9.61501846017335im

Parameters:
n = 0.01 + 0.0im
λ = 1.031036435632696e-8 + 5.386760588610906e-10im
τ = 0.001 + 0.0im
Y = 1.0e15 + 0.0im
, 10, 1)
```

## Maximum iterations

As per default five iterations are performed. This iteration number can be
customized using the keyword `maxiter` For instance the following command
runs 30 iterations

```julia
sol,nn,flag = mslp(L, z, maxiter=30)
```

```
(####Solution####
eigval:
ω = 1710.6977772393425 + 9.615018460173417im

Parameters:
n = 0.01 + 0.0im
λ = 1.1058932625637208e-9 - 2.9019411781324465e-11im
τ = 0.001 + 0.0im
Y = 1.0e15 + 0.0im
, 30, 1)
```

# A stopping crtierion

Usually, the procedure converges after a few iteration steps and further
iterations do not significantly improve the solution. For this reason the
keyword `tol` allows for setting a threshold, such that two consecutive
iterates for the eigenvalue fulfill `abs(ω_0-ω_1)`, the iteration is
terminated. This keyword defaults to `tol=0.0` and, thus, `maxiter` iterations
will be performed. However, it is highly recommended to set the parameter to
some positive value whenever possible.

```julia
sol,nn,flag = mslp(L, z, maxiter=30, tol=1E-10)
```

```
(####Solution####
eigval:
ω = 1710.6977772393418 + 9.615018460170836im

Parameters:
n = 0.01 + 0.0im
λ = 7.344557969393061e-8 + 2.3398225911910384e-9im
τ = 0.001 + 0.0im
Y = 1.0e15 + 0.0im
, 5, 0)
```

Note that the tolerance is also a bound for the error associated with the
computed eigenvalue. Tip: The 64-bit floating numbers have an accuracy of
`eps≈1E-16`. Hence, the threshold `tol` should not be chosen less than 16
orders of magnitude smaller than the expected eigenvalue. Indeed, `tol=1E-10`
is already close to the machine-precision we can expect for the Rijke-tube
model.
## Convergence checks

Technically, the iteration may stop because of slow progress in the iteration
rather than actual convergence. A simple indicator for actual convergence is
the auxiliary eigenvalue `λ`. The closer it is to `0` the better is the
quality of the computed solution. To help the identification of falsely
terminated iterations you can specify another tolerance `lam_tol`. If
`abs(λ)>lam_tol` the computed eigenvalue is deemed to be spurious. This
feature only makes sense when a termination threshold has been specified.
As in the following example:

```julia
sol,nn,flag = mslp(L, z, maxiter=30, tol=1E-10, lam_tol=1E-8)
```

```
(####Solution####
eigval:
ω = 1710.6977772393418 + 9.615018460170836im

Parameters:
n = 0.01 + 0.0im
λ = 7.344557969393061e-8 + 2.3398225911910384e-9im
τ = 0.001 + 0.0im
Y = 1.0e15 + 0.0im
, 5, 2)
```

## Faster convergence
In order to improve the convergence rate, higher-order derivatives may be
computed. You high-order perturbation theory to improve the
iteration via the `order` keyword. For instance, third-order theory is used
in this example

```julia
sol,nn,flag = mslp(L, z, maxiter=30, tol=1E-10, lam_tol=1E-8, order=3)
```

```
(####Solution####
eigval:
ω = 1710.6977772393466 + 9.615018460173351im

Parameters:
n = 0.01 + 0.0im
λ = -1.7829967501986358e-8 + 2.73112255541197e-8im
τ = 0.001 + 0.0im
Y = 1.0e15 + 0.0im
, 3, 2)
```

Under the hood the routine calls ARPACK to utilize Arnoldi's method to solve
the linear (auxiliary) eigenvalue problem. This method can compute more than
one eigenvalue close to some initial guess. Per default, only one is sought
but via the `nev` keyword the number can be increased. This results
in an increased computation time but provides more candidates for the next
iteration, potentially improving the convergence.

```julia
sol,nn,flag = mslp(L, z, maxiter=30, tol=1E-10, lam_tol=1E-8, nev=3)
```

```
(####Solution####
eigval:
ω = 1710.6977772393418 + 9.615018460170836im

Parameters:
n = 0.01 + 0.0im
λ = 7.344557969393061e-8 + 2.3398225911910384e-9im
τ = 0.001 + 0.0im
Y = 1.0e15 + 0.0im
, 5, 2)
```

## Summary

Iterative solver like the `mslp` solver
are a great opportunity to refine solutions found from Beyn's integration-based
method. The tutorial gave you more insight in how to use the keywords for `mslp`.

The next tutorial will explain how to post-process highly accurate solutions
from an iterative solver using perturbation theory.

# References
[1] S. Güttel and F. Tisseur, The Nonlinear Eigenvalue Problem, 2017, <http://eprints.ma.man.ac.uk/2538/>

[2] G.A. Mensah, A. Orchini, J.P. Moeck, Perturbation theory of nonlinear, non-self-adjoint eigenvalue problems: Simple eigenvalues, JSV, 2020, [doi:10.1016/j.jsv.2020.115200](https://doi.org/10.1016/j.jsv.2020.115200)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

