```@meta
EditURL = "<unknown>/tutorial_04_perturbation_theory.jl"
```

# Tutorial 04 Perturbation Theory

## Introduction

This tutorial demonstrates how high-order perturbation theory is utilized using
the WavesAndEigenvalues package. You may ask: 'What is perturbation theory?'
Well, perturbation theory deals with the approximation of a solution of a
mathematical problem based on a *known* solution of a problem similar to the
problem of interest. (Puhh, that was a long sentence...) The problem is said
to be perturbed from a baseline solution.
You can find a comprehensive presentation of the internal algorithms in
[1,2].

!!! note
    The following example uses perturbation theory up to 30th order. Per default
    WavesAndEigenvalues.jl is build to run perturbation theory to 16th order. You
    can reconfigure your installation for higher orders by setting an environment
    variable and then rebuild the package. For instance this tutorial requires
    `ENV["JULIA_WAE_PERT_ORDER"]=30` for setting the variable and a subsequent
    `] build WavesAndEigenvalues` to rebuild the package. This process may take a
    few minutes. If you don't make the environment variable permanent,
    these steps will be necessary once every time you got any update to the
    package from Julia's package manager. Also note that higher orders take
    significantly more memory of your installation directory.
## Model set-up and baseline-solution
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
P0=101325.0 # ambient pressure in Pa
A=pi*0.025^2 # cross sectional area of the tube
Q02U0=P0*(Tb/Tu-1)*A*γ/(γ-1) #the ratio of mean heat release to mean velocity Q02U0
x_ref=[0.0; 0.0; -0.00101] #reference point
n_ref=[0.0; 0.0; 1.00] #directional unit vector of reference velocity
n=1 #interaction index
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
n	1.0 + 0.0im
λ	Inf + 0.0im
ω	0.0 + 0.0im
τ	0.001 + 0.0im
Y	1.0e15 + 0.0im

```

To obtain a baseline solution we solve it using an iterative solver with
a very small stopping tollerance of `tol=1E-11`.

```julia
sol,nn,flag=mslp(L,340*2*pi,maxiter=20,tol=1E-11)
```

```
(####Solution####
eigval:
ω = 1075.325211506839 + 372.1017670372039im

Parameters:
n = 1.0 + 0.0im
λ = 2.336203477100084e-8 + 6.515535987345233e-10im
τ = 0.001 + 0.0im
Y = 1.0e15 + 0.0im
, 8, 0)
```

this small tolerance is necessary because as the name suggests the base-line
solution is the basis to our subsequent steps. If it is inaccurate we will
definetly also encounter inaccurate approximations to other configurations
than the base-line set-up.
## Taylor series
Probably, you are somewhat familiar with the concept of Taylor series
which is a basic example of perturbation theory. There is some function
``f`` from which the value ``f(x_0)`` is known together with some derivatives at the
same point, the first ``N`` say. Then, we can approximate ``f(x_0+\Delta)`` as:
```math
f(x_0+\Delta)\approx \sum_{n=0}^N f_n(x_0)Δ^n
```

Let's consider the time delay `τ`. Our problem depends exponentially on this
value.  We can utilize a fast perturbation algorithm to compute the first
20 Taylor-series coefficients of the eigenfrequency `ω` w.r.t. `τ` by just
typing

```julia
perturb_fast!(sol,L,:τ,20)
```

The output shows you how long it takes to compute the coefficients.
Obviously, it takes longer and longer for higher order coefficients.
Note, that there is also an algorithm `perturb!(sol,L,:τ,20)` which does
exactly the same as the fast algorithm but is slower, when it comes to high
orders (N>5).

Both algorithms populate a field in the solution object holding the Taylor
coefficients.

```julia
sol.eigval_pert
```

```
Dict{Symbol,Any} with 1 entry:
  Symbol("τ/Taylor") => Complex{Float64}[1075.33+372.102im, -2.62868e5+3.40796e5im, -1.79944e8-1.475e8im, 9.4741e10-1.57309e11im, 1.66943e14+8.14274e13im, -8.3483e16+1.86246e17im, -2.15622e20-9.17653e19im, 1.05357e23-2.58704e23im, 3.19354e26+1.25588e26im, -1.54315e29+4.02817e29im, -5.16822e32-1.94111e32im, 2.48748e35-6.72431e35im, 8.85193e38+3.23647e38im, -4.26474e41+1.17688e42im, -1.57803e45-5.68031e44im, 7.63541e47-2.13155e48im, 2.89779e51+1.0345e51im, -1.41133e54+3.9619e54im, -5.44414e57-1.93716e57im, 2.67322e60-7.51468e60im, 1.04148e64+3.70668e63im]
```

We can use these values to form the Taylor-series approximation and for
convenience we can just do so by calling to the solution object itself.
For instance let's assume we would like to approximate the value of the
eigenfrequency when `τ==0.0015` based on the 20th order series expansion.

```julia
Δ=0.0005
ω_approx=sol(:τ,τ+Δ,20)
```

```
916.7085040155473 + 494.3258317478708im
```

Let's compare this value to the true solution:

```julia
L.params[:τ]=τ+Δ #change parameter
sol_exact,nn,flag=mslp(L,ω_approx,maxiter=20,tol=1E-11) #solve again
ω_exact=sol_exact.params[:ω]
println(" exact=$(ω_exact/2/pi)  vs  approx=$(ω_approx/2/pi))")
```

```
Launching MSLP solver...
Iter   dz:     z:
----------------------------------
0			Inf	916.7085040155473 + 494.3258317478708im
1			0.006009923194378605	916.7036137652204 + 494.32932526022125im
2			2.5654097406043415e-8	916.7036137579167 + 494.3293252848137im
3			1.622588396708824e-11	916.7036137579236 + 494.329325284799im
4		5.804108531973829e-9	2.157056083093525e-12	916.7036137579256 + 494.32932528479967im
Solution has converged!
...finished MSLP!
#####################
 Results 
#####################
Number of steps: 4
Last step parameter variation:2.157056083093525e-12
Auxiliary eigenvalue λ residual (rhs):5.804108531973829e-9
Eigenvalue:916.7036137579256 + 494.32932528479967im
 exact=145.89791147977746 + 78.67495563435732im  vs  approx=145.89868978845095 + 78.6743996206862im)

```

Clearly, the approximation matches the first 4 digits after the point!

We can also compute the approximation at any lower order than 20.
For instance, the first-order approximation is:

```julia
ω_approx=sol(:τ,τ+Δ,1)
println(" first-order approx=$(ω_approx/2/pi)")
```

```
 first-order approx=150.22496667319837 + 86.34150633981955im

```

Note that the accuracy is less than at twentieth order.

Also note that we cannot compute the perturbation at higher order than 20
because we have only computed the Taylor series coefficients up to 20th order.
Of course we could, prepare higher order coefficients, let's say up to 30th
order by

```julia
perturb_fast!(sol,L,:τ,30)
```

and then

```julia
ω_approx=sol(:τ,τ+Δ,30)
println(" 30th-order approx=$(ω_approx/2/pi)")
```

```
 30th-order approx=145.8978874014616 + 78.67497208762059im

```

Ok, we've seen how we can compute Taylor-series coefficients and how to
evaluate them as a Taylor series. But how do we choose parameters like
the baseline set-up or the perturbation order to get a reasonable estimate?
Well there is a lot you can do but, unfortunately, there are a lot of
misunderstandings when it comes to the quality perturbative approximations.

Take a second and try to answer the following questions:
1. What is the range of Δ in which you can expect reliable results from the approximation, i.e., the difference from the true solution is small?
2. How does the quality of your approximation improve if you increase the number of known derivatives?
3. What else can you do to improve your solution?

Regarding the first question, many people misbelieve that the approximation is
good as long as Δ (the shift from the expansion point) is small. This is right
and wrong at the same time. Indeed, to classify Δ as small, you need to
compare it against some value. For Taylor series this is the radius of
convergence. Convergence radii can be quite huge and sometimes super small.
Also keep in mind that the numerical value of Δ in most engineering
applications is meaningless if it is not linked to some unit of measure.
(You know the joke: What is larger, 1 km or a million mm?)

So how do we get the radius of convergence? It's as simple as

```julia
r=conv_radius(sol,:τ)
```

```
30-element Array{Float64,1}:
 0.0026438071359421856
 0.0018498027477203886
 0.0012670310435927681
 0.0009886548876531008
 0.0009100554815832927
 0.0008709697017278116
 0.000838910762099998
 0.0008140058757174155
 0.0007955250149644752
 0.0007813536279922974
 0.0007700125089152769
 0.0007607027504841114
 0.000752936147548007
 0.0007463656620716248
 0.0007407359084332154
 0.0007358585624584575
 0.0007315926734366027
 0.0007278304643904657
 0.0007244879957019495
 0.0007214989316859786
 0.0007188101801265639
 0.0007163787543387638
 0.0007141694816882979
 0.0007121533078940557
 0.0007103060243302641
 0.0007086072999209774
 0.000707039935855723
 0.0007055892855979911
 0.0007042427990056821
 0.0007029896606802446
```

The return value here is an array that holds N-1 values, where N is the order
of Taylor-series coefficients that have been computed for `:τ`. This is
because the convergence radius is estimated based on the ratio of two
consecutive Taylor-series coefficients. Usually, the last entry of r should
give you the best approximate.

```julia
println("Best estimate available: r=$(r[end])")
```

```
Best estimate available: r=0.0007029896606802446

```

However, there are cases where the estimation procedure for the convergence
radius is not appropriate. That is why you can have a look on the other entries
in r to judge whether there is smooth convergence.
Let's see what's happening if we evaluate the Taylor series beyond it's radius
of convergence. We therefore try to find an approximation for a two that is
 twice our estimate for r away from the baseline value.

```julia
ω_approx=sol(:τ,τ+r[end]*2,20)
L.params[:τ]=τ+r[end]*2
sol_exact,nn,flag=mslp(L,ω_approx,maxiter=20,tol=1E-11) #solve again
ω_exact=sol_exact.params[:ω]
println(" exact=$(ω_exact/2/pi)  vs  approx=$(ω_approx/2/pi))")
```

```
Launching MSLP solver...
Iter   dz:     z:
----------------------------------
0			Inf	8.879592083136003e6 - 1.1737122716657885e6im
1			4.0600651883038566e6	4.870901607915898e6 - 529872.223861651im
2			1.682178480929161e6	3.225729888267407e6 - 178967.03072543198im
3			437944.1873905101	2.8166148076510336e6 - 22698.153205330833im
4			34045.990823039225	2.7910564230047897e6 - 205.96644198751892im
5			207.55279886998233	2.7910300114278947e6 - 0.1009691061235003im
6			0.012185956885981062	2.791030020923184e6 - 0.1086069737550432im
7			4.672416133789498e-10	2.7910300209231847e6 - 0.10860697371664671im
8		0.0011653820248254138	2.2028212587343887e-13	2.7910300209231847e6 - 0.10860697371686699im
Solution has converged!
...finished MSLP!
#####################
 Results 
#####################
Number of steps: 8
Last step parameter variation:2.2028212587343887e-13
Auxiliary eigenvalue λ residual (rhs):0.0011653820248254138
Eigenvalue:2.7910300209231847e6 - 0.10860697371686699im
 exact=444206.22414780094 - 0.01728533672129094im  vs  approx=1.413230972670755e6 - 186802.10980322777im)

```

Outch, the approximation is completely off. Also note that the exact solutions
has a ridicoulously high value. What's happening here? Well, because  we are
evaluating the power series outside of its convergence radius, we find an
extremely high value. Keep in mind that any quadratic or higher polynomial will
tend to infinity when its argument tends to infinity. Because our estimate is
nonsense but high and because we use this estimate to initialize our iterative
solver, we also find a high exact eigenvalue which of course does not match our
estimate.
Remember question 2? Maybe the estimate gets better if we increase the
perturbation order ? At least the last time we've seen an improvement.
But this time...

```julia
ω_approx=sol(:τ,τ+r[end]+0.001,30)
println("approx=$(ω_approx/2/pi))")
```

```
approx=-4.308053889489341e11 + 1.7221050696555923e10im)

```

things get *way* worse! A higher order polynomial tends faster to infinity,
therefore our is now extremely large. Indeed, if we would feed it into the
iterative solver it would trigger an error, because the code can't handle numbers
this large.  The large magnitude of the result is exactly the reason why some
clever person named r the radius of *convergence*. Beyond that radius we cannot
make the Taylor-series converge without shifting the expansion point. You might
think: "*Well, then let's shift the expansion point!*" While this is definitely
a valid approach, there is something better you can do...
## Series accelartion by Padé approximation

The Taylor-series cannot converge beyond its radius of convergence. So why
not trying something else than a Taylor-series? The truncated Taylor series is
a polynomial approximation to our unknown relation ω=ω(τ). An alternative
(if not to say a generalization) of this approach is to consider rational
polynomial approximations.  So instead of
$f(x_0+Δ)≈∑_{n=0}^N f_n(x_0)Δ^n$
we try something like
```math
f(x_0+Δ)≈\frac{\sum_{l=0}^L a_l(x_0)Δ^l }{ 1 + \sum_{m=0}^M b_m(x_0)Δ^m }
```
This approach is known as Padé approximation.
In order to get the same asymptotic behavior close to our expansion point, we
demand that the truncated Taylor series expansion of our Padé approximant is
identical to the approximation of our unknown function. Here comes the clou:
we already know these values because we computed the Taylor-series
coefficients. Only little algebra is needed to convert the Taylor coefficients
to Padé coefficients and all of this is build in to the solution type. All you
need to know is that the number of Padé coefficients is related to the number
of Taylor coefficients by the simple formula L+M=N. Let's give it a try and
compute the Padé approximant for L=M=10 (a so-called diagonal Padé
approximant). This is just achieve by an extra argument for the call to our
solution object.

```julia
ω_approx=sol(:τ,τ+r[end]*2,10,10)
sol_exact,nn,flag=mslp(L,ω_approx,maxiter=20,tol=1E-11)
ω_exact=sol_exact.params[:ω]
println(" exact=$(ω_exact/2/pi)  vs  approx=$(ω_approx/2/pi))")
```

```
Launching MSLP solver...
Iter   dz:     z:
----------------------------------
0			Inf	668.5373997804821 + 529.4636751544649im
1			4.96034674936746e-7	668.537399929806 + 529.4636746814398im
2		1.1965642295134598e-8	4.1740258326821776e-12	668.537399929804 + 529.4636746814361im
Solution has converged!
...finished MSLP!
#####################
 Results 
#####################
Number of steps: 2
Last step parameter variation:4.1740258326821776e-12
Auxiliary eigenvalue λ residual (rhs):1.1965642295134598e-8
Eigenvalue:668.537399929804 + 529.4636746814361im
 exact=106.40103184063163 + 84.26676101314976im  vs  approx=106.40103181686632 + 84.26676108843462im)

```

Wow! There is now only a difference of about 0.00001hz between the true solution
and it's estimate. I'would say that's fine for many engineering applications.
What's your opinion?
## Summary
You learnt how to use perturbation theory. The main computational costs for
this approach is the computation of Taylor-series coefficients. Once you have
them everything comes down to ta few evaluations of scalar polynomials.
Especially, when you like to rapidly evalaute your model for a bunch of values
in a given parameter range. This approach should speed up your computations.

The next tutorial will familiarize you with the use of higher order
finite elements.

# References
[1] G.A. Mensah, Efficient Computation of Thermoacoustic Modes, Ph.D. Thesis, TU Berlin, 2019, [doi:10.14279/depositonce-8952 ](http://dx.doi.org/10.14279/depositonce-8952)

[2] G.A. Mensah, A. Orchini, J.P. Moeck, Perturbation theory of nonlinear, non-self-adjoint eigenvalue problems: Simple eigenvalues, JSV, 2020, [doi:10.1016/j.jsv.2020.115200](https://doi.org/10.1016/j.jsv.2020.115200)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

