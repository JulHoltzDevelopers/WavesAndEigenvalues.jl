# # Exercise 1: Nonlinear Eigenvalue Solver
#
# This exercise aims at introdcing the package *WavesAndEigenvalues* applied to the
# identification of thermoacoustic instabilities. You will learn how to import and use a mesh,
# build a parametric dependent thermoacoustic operator, and use adjoint methods to solve
# the resulting nonlinear eigenvalue problems.
#
# ## Introduction
#
# The Helmholtz equation is the Fourier transform of the acoustic wave equation.
# With the presence of an unsteady heat release source, it reads:
#
# $$∇⋅(c² ∇ p̂) + ω² p̂ = -iω(γ-1)q̂$$
#
# The conventioned used to map frequency and time domains is $f'(t) --> f̂(ω)exp(+iωt)$.
# http://www.youtube.com/watch?v=g0NXlnsfqt0&t=43m6s
#
# Here $c$ is the speed-of-sound field, $p̂$ the (complex) pressure fluctuation modeshape,
# $ω$ the (complex) frequency of the problem, $i$ the imaginary unit, $γ$ the ratio of
# specifc heats, and $q̂$ the (frequency dependent) heat release fluctuation response.
#
# Together with some boundary conditions the Helmholtz equation models
# thermoacoustic instability in a cavity. The minimum required inputs to specify
# the problem are
#
# For the acoustics:
# 1. the 3D shape of the domain;
# 2. the speed-of-sound field in the domain
# 3. the boundary conditions
#
# For an (active) flame, q̂≠0:
# 4. some gas porperties;
# 5. and a flame response model that links q̂ to p̂
#
## #jl
#
# ## Finite Element discretization of the Helmholtz equation
#
# To start, we shall consider a straight duct of length $L = 0.5$~m and radius $r=0.025$~m,
# filled with air.  The gas properties upstream of the flame are $\gamma=1.4$, $\rho=1.17$,
# and $p_0=101325$, and $T_u=300$. We will set the boundary conditions to be closed (upstream)
# and opened (downstream). For this simple system, you may know and/or be able to calculate the
# acoustic eigenvalues analytically. This fundamental frequency is $f_1 = c/4L$. We will benchmark
# the FE solver results with these known values.
#
# To solve this problem with WavesAndEigenvalues you need a mesh, the equations, a discretization scheme, and an eigenvalue solver. The mesh, called "Rijke{\textunderscore}mm.msh" is  is given to you. All the rest is pre-coded for you in WavesAndEigenvalues, you only need to define the problem.
#
# To learn how to load a mesh, define the equations to be solved, and build a discretization matrix, use the documentation of the package.
#
#
# ### The Helmholtz equation
using WavesAndEigenvalues.Helmholtz
#
# To define the problem, we need to specify what equations are to be solved in each
# part of the domain. In a finite element formulation, these are expressed in terms of
# a discretized version of the equations that govern the linear thermoacoustic dynamics
# expressed in *weak form*.
#
# WavesAndEigenvalues contains pre-defined standard finite element discretization matrices.
# These include e.g., mass and stiffness matrices, boundary terms, etc.
# You can find these matrices in the src/FEM.jl file.
#
# We shall start by considering a Rijke tube, and then move to more complex geometries.
# Open the pdf file "Exercise1".
#
# First you will need to load the Helmholtz solver. The following line brings all
# necessary tools into scope:
#
#
#
# A Rijke tube is the most simple thermo-acoustic configuration. It comprises a
# tube with an unsteady heat source inside the tube. This excercise will
# walk you through the basic steps of setting up a Rijke tube in a
# Helmholtz-solver stability analysis.




## #jl
# ## Mesh
# Now we can load the mesh file. It is the essential piece of information
# defining the domain shape . This tutorial's mesh has been specified in mm
# which is why we scale it by `0.001` to load the coordinates as m:
mesh=Mesh("Rijke_mm.msh",scale=0.001)

# You can have a look at some basic mesh data by printing it
print(mesh)
# In an interactive session, is suffices to just type the mesh's variable name
# without the enclosing `print` command.
#
#
# This info tells us that the mesh features `1006` points as vertices
# `1562` (addressable) triangles on its surface and `3380`tetrahedra forming
# the tube. Note, that there are currently no lines stored. This is because
# line information is only needed when using second-order finite elements. To
# save memory this information is only assembled and stored when explicitly
# requested. We will come back to this aspect in a later tutorial.
#
# You can also see that the mesh features several named domains such as
# `"Interior"`,`"Flame"`, `"Inlet"` and `"Outlet"`. These domains are used to
# specify certain conditions for our model.
#
# A model descriptor is always given as a dictionairy featureing domain names
# as keys and tuples as values. The first tuple entry is a Julia symbol
# specifying the operator to be build on the associated domain, while the second
# entry is again a tuple holding information that is specific to the chosen
# operator.

## #src
# ## Model set-up
# For instance, we certainly want to build the wave operator on the entire mesh.
# So first, we initialize an empty dictionairy
dscrp=Dict()

# And then specify that the wave operator should be build everywhere. For the
# given mesh the domain-specifier to address the entire resononant cavity is
# `"Interior"` and in general the symbol to create the wave operator is
# `:interior`. It requires no options so the options tuple remains empty.
dscrp["Interior"]=(:interior, ())

## src
# ## Boundary Conditions
# The tube should have a pressure node at its outlet, in order to model an
# open-end boundary condition. Boundary conditions are specified in terms of
# admittance values. A pressure note corresponds to an infinite admittance.
# To represent this on a computer we just give it a crazily high value like
# `1E15`. We will also need to specify a variable name for our admmittance value.
# This will allow to quickly change the value after discretization of the
# model by addressing by this very name. This feature is one of the core
# concepts of the solver. As will be demonstrated later.
#
# The complete description of our boundary condition reads
dscrp["Outlet"]=(:admittance, (:Y,1E15))

# You may wonder why we do not specify conditions for the other boundaries of
# the model. This is because the discretization is based on Bubnov-Galerkin
# finite elements. All unspecified, boundaries will therefore *naturally* be
# discretized as solid walls.

## #src
# ## Flame
# The Rijke tube's main feauture is a domain of heat release. In the current
# mesh there is a designated domain `"Flame"` addressing a thin volume at the
# center of the tube. We can use this key to specify a flame with simple
# n-tau-dynamics. Therefore, we first need to define some basic gas properties.

γ=1.4 #ratio of specific heats
ρ=1.17 #density at the reference location upstream to the flame in kg/m^3
Tu=300.0    #K unburnt gas temperature
Tb=1200.0    #K burnt gas temperature
P0=101325.0 # ambient pressure in Pa
A=pi*0.025^2 # cross sectional area of the tube
Q02U0=P0*(Tb/Tu-1)*A*γ/(γ-1) #the ratio of mean heat release to mean velocity Q02U0

# We will also need to provide the position where the reference velocity has
# been taken and the direction of that velocity
x_ref=[0.0; 0.0; -0.00101]
n_ref=[0.0; 0.0; 1.00] #must be a unit vector
# And of course, we need some values for n and tau
n=0.0 #interaction index
τ=0.001 #time delay
# With these values the specification of the flame reads
dscrp["Flame"]=(:flame,(γ,ρ,Q02U0,x_ref,n_ref,:n,:τ,n,τ))
# Note that the order of the values in the options tuple is important. Also note
# that we assign the symbols `:n`and `:τ` for later analysis.

## #jl
# ## Speed of Sound
# The description of the Rijke tube model is nearly complete. We just need to
# specify the speed of sound field.  For this example, the field is fairly
# simple and can be specified analytically, using the `generate_field` function
# and a custom three-parameter function.
R=287.05 # J/(kg*K) specific gas constant (air)
function speedofsound(x,y,z)
    if z<0.
        return sqrt(γ*R*Tu)#m/s
    else
        return sqrt(γ*R*Tb)#m/s
    end
end
c=generate_field(mesh,speedofsound)
# Note that in more elaborate models, you may read the field `c` from a file
# containing simulation or experimental data, rather than specifying it
# analytically.

## #src
# ## Model Discretization
# Now we have all ingredients together and we can discretize the model.
L=discretize(mesh,dscrp,c)
# The return value `L` here is a family of linear operators. You can display
# an algebraic representation of the family plus a list of associated parameters
# by just printing it
print(L)
# (In an interactive session the enclosing print function is not necessary.)
#
# You might notice that the list contains all parameters that have been
# specified during the model description (`n`,`τ`,`Y`). However, there are also
# two parameters that were added by default: `ω` and `λ`. `ω` is the complex
# frequency of the model and `λ` an auxiliary value that is important for the
# eigenfrequency computation.
## #src
# ## Solving the Model
# ### Global Solver
# We can solve the model for some eigenvalues using one of the eigenvalue
# solvers. The package provides you with two types of eigenvalue solvers. A global
# contour-integration-based solver. That finds you all eigenvalues inside of a
# specified contour Γ in the complex plane.
Γ=[150.0+5.0im, 150.0-5.0im, 1000.0-5.0im, 1000.0+5.0im].*2*pi #corner points for the contour (in this case a rectangle)
Ω, P = beyn(L,Γ,l=10,N=64, output=true)

# The found eigenvalues are stored in the array `Ω`. The corresponding
# eigenvectors are the columns of `P`.
# The huge advantage of the global eigenvalue solver is that it finds you
# multiple eigenvalues. However, its accuracy may be low and sometimes it
# provides you with outright wrong solutions.
# However, for the current case the method works as we can varify that in our
# search window there are two eigenmodes oscilating at 272 and 695 Hz
# respectively:
for ω in Ω
    println("Frequency: $(ω/2/pi)")
end
## #src
# ### Local Solver
# To get high accuracy eigenvalues , there is also local iterative eigenvalue
# solver. Based on an initial guess, it only finds you one eigenvalue at a time
# but with machine-precise accuracy.
sol,nn,flag=householder(L,250*2*pi,output=true);

# The return values of this solver are a solution object `sol`, the number of
# iterations `nn` performed by the solver and an error flag `flag` providing
# information on the quality of the solution. If `flag>0` the solution has
# converged. If `flag<0`something went wrong. And if `flag==0`... well, probably
# the solution is fine but you should check it.
#
# In the current case the flag is -1 indicating that the maximum number of
# iterations has been reached. This is because we haven't specified a stopping
# criterion and the iteration just ran for the maximum number of iterations.
# Nevertheless, the 272 Hz mode is correct. Indeed, as you can see from the
# printed output it already converged to machine-precision after 5 iterations.

## #src
# ## Changing a specified parameter.
#
# Remember that we have specified the parameters `Y`,`n`, and `τ`?. We can change
# these easily without rediscretizing our model. For instance currently `n==0.0`
# so the solution is purely accoustic. The effect of the flame just enters the
# problem via the speed of sound field (this is known as passive flame). By
# setting the interaction index to `1.0` we activate the flame.
L.params[:n]=1
# Now we can just solve the model again to get the active solutions
sol_actv,nn,flag=householder(L,245*2*pi-82im*2*pi,output=true,order=1);
# The eigenfrequency is now complex:
sol_actv.params[:ω]/2/pi
# with a growth rate `-imag(sol_actv.params[:ω]/2/pi)≈  59.22`

# Instead of recomputing the eigenvalue sol_acby one of the solvers. We can also
# approximate the eigenvalue as a function of one of the model parameters
# utilizing high-order perturbation theory. For instance this gives you
# the 30th order diagonal Padé estimate expanded from the passive solution.

perturb_fast!(sol,L,:n,16) #compute the coefficients
freq(n)=sol(:n,n,8,8)/2/pi #create some function for convenience
freq(1) #evaluate the function at any value for `n` you are interested in.

# The method is slow when only computing one eigenvalue. But its computational
# costs are almost constant when computing multiple eigenvalues. Therefore,
# for larger parametric studies this is clearly the method of choice.

## #src
# ## VTK Output for Paraview
#
# To visualizte your solutions you can store them as `"*.vtu"` files to open
# them with paraview. Just, create a dictionairy, where your modes are stored as
# fields.
data=Dict()
data["speed_of_sound"]=c
data["abs"]=abs.(sol_actv.v)/maximum(abs.(sol_actv.v)) #normalize so that max=1
data["phase"]=angle.(sol_actv.v)

vtk_write("tutorial_01", mesh, data) # Write the solution to paraview

# The `vtk_write` function automatically sorts your data by its interpolation
# scheme. Data that is constant on a tetrahedron will be written to a file
# `*_const.vtu`, data that is linearly interpolated on a tetrahedron will go
# to `*_lin.vtu`, and data that is quadratically interpolated to `*_quad.vtu`.
# The current example uses constant speed of sound on a tetrahedron and linear
# finite elements for the discretization of p. Therefore, two files are
# generated, namely `tutorial_01_const.vtu` containing the speed-of-sound field
# and `tutorial_01_lin.vtu` containing the mode shape. Open them with paraview
# and have a look!.

## src
# ## Summary
#
# You learnt how to set-up a simple Rijke-tube study. This  already introduced
# all basic steps that are needed to work with the Helmholtz solver. However,
# there are a lot of details you can fine tune and even features that weren't
# mentioned in this tutorial. The next tutorials will introduce these aspects in
# more detail.
