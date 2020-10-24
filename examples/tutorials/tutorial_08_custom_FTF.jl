# The model is the same Rijke tube configuration as in Tutorial 01. However,
# this time the FTF is customly defined. Let's start with the definition of the
# flame transfer function. The custom function must return the FTF and all its
# derivatives w.r.t. the complex freqeuncy. Therefore, the interface **must** be
# a function with two inputs: a complex number and an integer. For an n-τ-model
# it reads.
function FTF(ω::ComplexF64, k::Int=0)
    return n*(-1.0im*τ)^k*exp(-1.0im*ω*τ)
end
# Note, that I haven't defined `n` and `τ` here, but Julia is a quite generous
# programming language; it won't complain if we define them before making our
# first call to the function.

# Of course, it is sometimes hard to find a closed form expression for an FTF
# and all its derivatives. You need to specify at least these derivative orders
# that will be used by the methods applied to analyse the model. This is at
# least the first order derivative for beyn and the standard householder method.
# You might return derivative orders up to the order you need it using an if
# clause for each order. Anyway, the function may get complicated. The implicit
# method explained below is then a viable alternative.

# Now let's set-up the Rijke tube model. The data is the same as in tutorial 1
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
R=287.05 # J/(kg*K) specific gas constant (air)
speedofsound(x,y,z) = z<0. ? sqrt(γ*R*Tu) : sqrt(γ*R*Tb)
c=generate_field(mesh,speedofsound)
# btw this is the moment where we specify the values of n and τ
n=0.01 #interaction index
τ=0.001 #time delay

##
# but instead of passing n and τ and its values to the flame description we
# just pass the FTF
dscrp["Flame"]=(:flame,(γ,ρ,Q02U0,x_ref,n_ref,FTF)) #flame dynamics

#... discretize ...
L=discretize(mesh,dscrp,c)
# and done!
## Check the model
print(L)
# Independently of how you name it, your FTF will be displayed as `FTF`
# in the signature of the operator.

##
# Solving the problem shows that we get the same result as with the built-in
# n-τ-model in tutorial 1
sol,nn,flag=householder(L,340*2*pi,maxiter=20,tol=1E-11)
## changing the parameters
# Because the system parameters n and τ are defined in this outer scope,
# changing them will change the FTF and thus the model without rediscretization.
# For instance, you can completely deactivate the FTF by setting n to 0.
n=0
# Indeed the model now converges to purely acoustic mode
sol,nn,flag=householder(L,340*2*pi,maxiter=20,tol=1E-11)

# However, be aware that there is currently no mechanism keeping track of these
# parameters. It is the programmer's (that's you) responsibility to know the
# specification of the FTF when `sol` was computed.

# Also, note the parameters are defined in this scope. You cannot change it by
# iterating it in a for-loop or from within a function, due to julia's scoping
# rules.

## Something else then n-τ...
# Of course this is an academic example. In practice you will specify
# something very custom as FTF. Maybe a superposition of σ-τ-models or even
# a statespace model fitted to experimental data.

#TODO: implicit FTF modelling, vectorfit
