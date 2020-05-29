# #Tutorial 02
#
# ##Introduction
#
# In Tutorial 01 you learnt the basics of the Helmholtz solver. This tutorial
# will make you more familiar with the global eigenvalue solver. The solver
# implements an algorithm invented by W.-J. Beyn in [1]. The implementation
# closely follows [2].
# ### References
# [1] W.-J. Beyn, An integral method for solving nonlinear eigenvalue problems, Linear Algebra and its Applications, 2012, 436(10), p.3839-3863, https://doi.org/10.1016/j.laa.2011.03.030
# [2] P.E. Buschmann, G.A. Mensah, J.P. Moeck, Solution of Thermoacoustic Eigenvalue Problems with a Non-Iterative Method, J. Eng. Gas Turbines Power, Mar 2020, 142(3): 031022 (11 pages) https://doi.org/10.1115/1.4045076

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
## #jl
# The global eigenvalue solver has a range of parameters. Mandatory, are only
# the linear operator family `L` you like to solve and the contour `Γ`.
#
# The contour is specified as a list of complex numbers. Each of which defining
# a vertex of the polygon you want to scan for eigenvalues.
# In Tutorial 01 we specified the contour as rectangle by
Γ=[150.0+5.0im, 150.0-5.0im, 1000.0-5.0im, 1000.0+5.0im].*2*pi
# The verteces are traversed in the order you specify them. However, the
# direction does not matter, meaning that equivalently you can specify `Γ` as:
Γ=[1000.0+5.0im, 1000.0-5.0im, 150.0-5.0im, 150.0+5.0im]
# Note that you should *not* specify the first point as the last point.
# This will be automatically done for you by the solver.
#
# Again, Γ is a polygon, so, e.g., specifying only three points will make the
# search region a triangle. And to get curved regions you just approximate them
# by a lot of vertices. For instance, this polygon is a good approximation of
# a circle.
radius=10
center_point=100
Γ_circle=radius.*[cos(α)+sin(α)*1.0im for α=0:2*pi/1000:(2*pi-2*pi/1000)].+center_point
#
# Note that your contour does not need to be convex!
#
# ## Number of eigenvalues
#
# Beyn's algorithm needs an initial guess on how many eigenvalues `l` you expect
# inside your contour. Per default this is `l=5`. If this guess is less than
# the actual eigenvalues. The algorithms will miserably fail to compute any
# eigenvalue inside your contour correctly. You may, therefore, choose a large
# number for `l` to be on the save site. However, this will dramatically
# increase your computation time.
#
# There are several strategies to deal with this problem and each of which might
# be well suited in a certain situation:
# 1. Split your domain in several smaller domains. This reduces the potential
# number of eigenvalues in each of the subdomains.
# 2. Just repeatedly, rerun the solver with increasing values of `l` until
# the found eigenvalues do not change anymore.
#
# ##Qudrature points
#
# As Beyn's algorithm is based on contour integration a numerical quadrature
# method is performed under the hood. More precisely, the code performs a
# Gauss-Legendre integration on each of the edges of the Polygon `Γ`. The number
# of quadrature points for each of this integrations is `N=32` per default. But
# you can specify it as an optional parameter when calling the solver.
# For instance as the circular contour is already by by 1000 points. It is not
# necessary to utilize `N=16` quadrature points on each of the edges. Let's say
# `N=4` should be fine. Then you would call the solver like
Ω,P=beyn(L,Γ_circle,N=4,output=true)
## #src
# ## test singular values
#
# One step in Beyn's algorithm requires a singular value decomposition. Non-zero
# singular values and singular vectors associated with these are meant to be
# disregarded. Unfornately, on a computer zero could also mean a very small
# number, and *small* here is problem dependend. The code will disregard any
# singular value `σ<tol` where `tol = 0.0` Per default. This means, that per
# default only singular values that are exactly 0 will be disregraded. It is
# very unlikely that this will be the case for any singular value. You may
# specify you own threshold by the optional key-word argument `tol`
Ω,P=beyn(L,Γ, tol=1E-10)
# but this is rather an option for experts. Even the authors of the code do not
# know a systematic way of well defining this threshold for given configuration.
## #src
# ##Test the position of the eigenvalues
#
# Because there could be a lot of spurious eigenvalues as there is no correct
# implementation of the singular value test. A very simple test to chek whether
# your eigenvalues are correct is to see whether they are inside your contour
# `Γ`. This test is performed per default, but you can disable it using the
# keyword `pos_test`.
Ω,P=beyn(L,Γ, pos_test=false)
# Note, that if you disable the position test and do not specify a threshold for
# the singular values you will allways be returned with `l` eigenvalues by
# Beyn's algorithm. (And sometimes even eigenvalues that are outside
# of your contour are fairly correct... )
##   #src
# ##Toggle the progressbar
#
# Especially when searching for many eigenvalues, in large domains the algorithm
# may take some time to finish. You can, therefore, optionally display a
# progress bar together with an estimate of how long it will take for the
# routine to finish. This is done using the optional keyword `output`.
Ω,P=beyn(L,Γ, output=true)

## TODO Caveats
#
# - contours close to eigenvalues
# - refine your solutions using householder
# - mention degeneracy
