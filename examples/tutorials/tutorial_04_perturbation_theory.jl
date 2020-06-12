# # Tutorial 04 Perturbation Theory
#
# ## Introduction
#
# This tutorial demonstrates how perturbation theory is utilized using the
# WavesAndEigenvalues package. You may ask: 'What is perturbation theory?'
# Well, perturbation theory deals with the approximation of a solution of a
# mathematical problem based on a *known* solution of a problem similar to the
# problem of interest. (Puhh, that was a long sentence...) The problem is said
# to be perturbed from a baseline solution.
# You can find a comprehensive presentation of the internal algorithms in
# [1,2,3].
## #jl
# ## Model set-up and baseline-solution
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

# To obtain a baseline solution we solve it using the housholder iteration with
# a very small stopping tollerance of `tol=1E-11`.
sol,nn,flag=householder(L,340*2*pi,maxiter=20,tol=1E-11)

# this small tolerance is necessary because as the name suggests the base-line
# solution is the basis to our subsequent steps. If it is inaccurate we will
# definetly also encounter inaccurate approximations to other configurations
# than the base-line set-up.
## #jl
# ## Taylor series
#
# Probably, you are somewhat familiar with the concept of Taylor series
# which is a basic example of perturbation theory. There is some function
# f from which the value f(x0) is known together with some derivatives at the
# same point, the first N say. Then, we can approximate f(x0+Δ) as:
#
# f(x0+Δ)≈∑_n=1^N f_n(x0)Δ^n

# Let's consider the time delay `τ`. Our problem depends exponentially on this
# value.  We can utilize a fast perturbation algorithm to compute the first
# 20 Taylor-series coefficients of the eigenfrequency `ω` w.r.t. `τ` by just
# typing
perturb_fast!(sol,L,:τ,20)
# The output shows you how long it takes to compute the coefficients.
# Obviously, it takes longer and longer for higher order coefficients.
# Note, that there is also an algorithm `perturb!(sol,L,:τ,20)` which does
# exactly the same as the fast algorithm but is slower, when it comes to high
# orders (N>5).

# Both algorithms populate a field in the solution object holding the Taylor
# coefficients.
sol.eigval_pert

# We can use these values to form the Taylor-series approximation and for
# convenience we can just do so by calling to the solution object itself.
# For instance let's assume we would like to approximate the value of the
# eigenfrequency when "τ==0.0" based on the 20th order series expansion.
ω_approx=sol(:τ,0.0,20)

# Let's compare this value to the true solution:
L.params[:τ]=0.0 #change parameter
sol_exact,nn,flag=householder(L,340*2*pi,maxiter=20,tol=1E-11) #solve again
ω_exact=sol_exact.params[:ω]
println(" exact=$(ω_exact/2/pi)  vs  approx=$(ω_approx/2/pi))")
# Clearly, the approximation matches the first 4 digits after the point!

# We can also compute the approximation at any lower order than 20.
# For instance, the first-order approximation is:
ω_approx=sol(:τ,0.0,1)
println(" first-order approx=$(ω_exact/2/pi)")

# Note that the accuracy is less than at twentieth order.
#
# Also note that we cannot compute the perturbation at higher order than 20
# because we have only computed the Taylor series coefficients up to 20th order.
# Of course we could, prepare higher order coefficients, let's say up to 30th
# order by
perturb_fast!(sol,L,:τ,30)
# and then
ω_approx=sol(:τ,0.0,30)
println(" 30th-order approx=$(ω_exact/2/pi)")

# Indeed, `30` is a special limit because per default the WavesAndEigenvalues
# package installs the perturbation algorithm on your machine up to 30th order.
# This is because higher orders would consume significantly more memory in
# the installation directory and also the computation of higher orders may take
# some time.  However, if you really want to use higher orders. You can rebuild
# your package by first setting a special environment variable, to your desired
# order -- say 40 -- by `ENV["JULIA_WAE_PERT_ORDER"]=40` and then run
# `Pkg.build("WavesAndEigenvalues")`
# Note, that if you don't make the environment variable permanent,
# this will be necessary once every time you got any update to the package
# from julias package manager.


#TODO: Padé, speed, conv_radius


# Unfornately, there are a lot of misunderstandings when it comes to the
# validity of such approximations. Take a second and try to answer the following
# questions:
# 1. What is the range of Δ in which you can expect reliable results from
# the approximation, i.e., the difference from the true solution is small?
# 2. How does the quality of your approximation improve if you increase the
# number of known derivatives?
# 3. What else can you do to improve your solution?
#
# Regarding the first question, many people misbelieve that the approximation is
# good as long as Δ is small. This is right and wrong at the same time. Indeed,
# to classify Δ as small, you need to compare it agains some value. For Taylor
# series this is the radius of convergence. Convergence radii can be quite huge
# and sometime super small. Also keep in mind that the nuemrical value of Δ in
# most engineering applications is meaningless if it is not linked to some unit
# of measure. (You know the joke: What is larger, 1 km or a million mm?)
