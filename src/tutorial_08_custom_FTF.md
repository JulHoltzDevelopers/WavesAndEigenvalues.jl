```@meta
EditURL = "<unknown>/tutorial_08_custom_FTF.jl"
```

# Tutorial 08 Custom FTF

The typical use case for thermoacoustic stability assessment will require
case-speicific data. e.g, measured impedance functions or flame transfer
functions. This tutorial explains how to specify custom functions in order to
include them in your model.

## Model set-up

The model is the same Rijke tube configuration as in Tutorial 01. However,
this time the FTF is customly defined. Let's start with the definition of the
flame transfer function. The custom function must return the FTF and all its
derivatives w.r.t. the complex freqeuncy. Therefore, the interface **must** be
a function with two inputs: a complex number and an integer. For an n-τ-model
it reads.

```@example tutorial_08_custom_FTF
function FTF(ω::ComplexF64, k::Int=0)::ComplexF64
    return n*(-1.0im*τ)^k*exp(-1.0im*ω*τ)
end
```

!!! note
    `n` and `τ` aren't defined yet, but Julia is a quite generous
    programming language; it won't complain if we define them before making our
    first call to the function.

Of course, it is sometimes hard to find a closed-form expression for an FTF
and all its derivatives. You need to specify at least these derivative orders
that will be used by the methods applied to analyse the model. This is at
least the first order derivative for beyn and the standard mslp method.
You might return derivative orders up to the order you need it using an if
clause for each order. Anyway, the function may get complicated. The implicit
method explained further below is then a viable alternative.

For convenience we can also specify a display signature, that is used to nicely
display our function when the discretized model is displayed in the REPL. Let's
say in our case we want the function to be named `"ntau(z)"` where `"z"` is a
placeholder for whatever the name of the argument will be. We achieve this by
adding a method to our FTF function that takes a symbol as input and returns a
string.

```@example tutorial_08_custom_FTF
function FTF(z::Symbol)::String
    return "ntau($z)"
end
```

Now let's set-up the Rijke tube model. The data is the same as in tutorial 1

```@example tutorial_08_custom_FTF
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
```

btw this is the moment where we specify the values of n and τ

```@example tutorial_08_custom_FTF
n=1 #interaction index
τ=0.001 #time delay

#
```

but instead of passing n and τ and its values to the flame description we
just pass the FTF

```@example tutorial_08_custom_FTF
dscrp["Flame"]=(:flame,(γ,ρ,Q02U0,x_ref,n_ref,FTF)) #flame dynamics

#... discretize ...
L=discretize(mesh,dscrp,c)
```

and done! Note how the function is currectly displayed as "natu(ω)"  if we
wouldn't have named it, the FTF will be displayed as `"FTF(ω)"`
in the signature of the operator.

## Checking the model

Solving the problem shows that we get the same result as with the built-in
n-τ-model in Tutorial 1

```@example tutorial_08_custom_FTF
sol,nn,flag=mslp(L,340,maxiter=20,tol=1E-9,scale=2*pi)
# changing the parameters
```

Because the system parameters n and τ are defined in this outer scope,
changing them will change the FTF and thus the model without rediscretization.
For instance, you can completely deactivate the FTF by setting n to 0.

```@example tutorial_08_custom_FTF
n=0
```

Indeed the model now converges to purely acoustic mode

```@example tutorial_08_custom_FTF
sol,nn,flag=mslp(L,340*2*pi,maxiter=20,tol=1E-11)
```

However, be aware that there is currently no mechanism keeping track of these
parameters. It is the programmer's (that's you) responsibility to know the
specification of the FTF when `sol` was computed.

Also, note the parameters are defined in this scope. You cannot change it by
iterating it in a for-loop or from within a function, due to julia's scoping
rules.

```@example tutorial_08_custom_FTF
# Something else than n-τ...
```

Of course this is an academic example. In practice you will specify
something very custom as FTF. Maybe a superposition of σ-τ-models or even
a statespace model fitted to experimental data.

A faily generic way of handling was intoduced in [1]. This approach is not
considering a specific FTF but parametrizes the problem in terms of the flame
response. Using perturbation theory it is than possible to get the
frequency as a function of the flame response. Closing the system with a
specific FTF is than possible, and the eigenfrequency is found by solving a
scalar equation! Here is a demonstration on how that works:

```@example tutorial_08_custom_FTF
dscrp["Flame"]=(:flame,(γ,ρ,Q02U0,x_ref,n_ref)) #no flame dynamics specified at all
H=discretize(mesh,dscrp,c)
η0=0
H.params[:FTF]=η0
sol_base,nn,flag=mslp(H,340*2*pi,maxiter=20,tol=1E-11)
```

The model H has no flame transfer function, but the flame response appears as a
parameter `:FTF`. We set this parameter to 0 and solved the model, i.e. we
computed a possive solution.  This solution will serve as a baseline. Note, that
it is not necessary to use a passive solution as baseline, a nonzero baseline
flame response would also work. Anyway, the baseline solution is used to build
a XX-order pade-approximation. For convenience we import some tools to handle
the algebra.

```@example tutorial_08_custom_FTF
#
using WavesAndEigenvalues.Helmholtz.NLEVP.Pade: Polynomial, Rational_Polynomial, prodrange, pade
#
perturb_fast!(sol_base,H,:FTF,32)
#
func=pade(Polynomial(sol_base.eigval_pert[Symbol("FTF/Taylor")]),15,15)
ω(z,k=0)=func(z-η0,k)
#
```

## Newton-Raphson for finding flame response
```math
\eta = FTF(\omega(\Delta\eta)}-\eta_0
```
From Newton-Raphson we find

Δη_{k+1}=Δη_{k+1}-\frac{FTF(ω(Δη_k)}-Δη_k}{FTF'(ω(Δη_k)}ω''(Δη_k)-Δη_k-1}

```@example tutorial_08_custom_FTF
#
n=1
τ=0.001
η=sol_base.params[:FTF]
for ii=1:10
    global omeg,η
    omeg=ω(η)
    η-=(FTF(omeg)-η)/(FTF(omeg,1)*ω(η,1)-1)
    println(ii,":",η)
end
println(ω(η)/2pi," η :",η," Res:",FTF(omeg)-η)
#
sol_exact,nn,flag=mslp(L,ω(η)/2pi,maxiter=20,tol=1E-11,output=false) #solve again
ω_exact=sol_exact.params[:ω]
println(" exact=$(ω_exact/2/pi)  vs  approx=$(ω(η)/2pi))")
#
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

