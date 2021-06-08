```@meta
EditURL = "<unknown>/tutorial_09_forcing.jl"
```

# Tutorial 09 Forcing

The tutorial demonstrates how to incorporate source terms into your Helmholtz
model. This allows for the computation of acoustic transfer functions and
similar qunatities. Currently, only boundary sources are supported.
This tutorial will use the outlet wall of the rijke tube mesh, as an acoustic
source. Such a case would occur when simulating acoustic excitations from
a loudspeaker membrane.


## Model set-up.

The model is the same Rijke tube configuration as in Tutorial 01.

The only difference is that we put a keyowrd `:speaker` together with some
parameters in our model description. The parameters are the excitation strength
and the speaker admittance.

In particular, the speaker forcing is derived from an impedance boundary condition
p̂ - (cZ)/(iω) * ∇p̂⋅n=A where p̂ is the pressure fluctuation amplitude, c the
local speed of sound, Z a normalized impedance, i the imaginary unit, ω an
angular frequency, n the outward pointing unit normal, and A the excitation
strength. This equation is equal to -iωY/c * p̂+∇p̂⋅n=-iωY/c*A, where the admittance
is Y=1/Z. Note that for the special case of a a sound-soft impedance (Z=0 aka
Y=∞) the excitation level A is identical to the pressure amplitude p̂ (p̂=A).

Ok let's see how it works. We start with the standard header known from tutorial 01.

```@example tutorial_09_forcing
using WavesAndEigenvalues.Helmholtz
mesh=Mesh("Rijke_mm.msh",scale=0.001) #load mesh
dscrp=Dict() #initialize model descriptor
dscrp["Interior"]=(:interior, ()) #define resonant cavity
```

Now comes the magic line that specifies the excitation level and speaker
impedance:

```@example tutorial_09_forcing
dscrp["Outlet"]=(:speaker, (:A,1,:Y,1E15)); #specify outlet BC
nothing #hide
```

The rest is the usual set-up...

```@example tutorial_09_forcing
γ=1.4; #ratio of specific heats
ρ=1.225; #density at the reference location upstream to the flame in kg/m^3
Tu=300.0;    #K unburnt gas temperature
Tb=1200.0;    #K burnt gas temperature
P0=101325.0; # ambient pressure in Pa
A=pi*0.025^2; # cross sectional area of the tube
Q02U0=P0*(Tb/Tu-1)*A*γ/(γ-1); #the ratio of mean heat release to mean velocity Q02U0
x_ref=[0.0; 0.0; -0.00101]; #reference point
n_ref=[0.0; 0.0; 1.00]; #directional unit vector of reference velocity
n=0.01; #interaction index
τ=0.001; #time delay
dscrp["Flame"]=(:flame,(γ,ρ,Q02U0,x_ref,n_ref,:n,:τ,n,τ)); #flame dynamics
R=287.05; # J/(kg*K) specific gas constant (air)
speedofsound(x,y,z) = z<0. ? sqrt(γ*R*Tu) : sqrt(γ*R*Tb);
c=generate_field(mesh,speedofsound);
nothing #hide
```

Note that when calling `discretize` the `source` keyword is set to `true` and
a second term `rhs` is returned. This is the discretization of the right hand
side.

```@example tutorial_09_forcing
L,rhs=discretize(mesh,dscrp,c,source=true,order=:quad);
nothing #hide
```

Like `L`, `rhs` is a linear operator family but note that it is a family of
vectors not matrices. You can reset parameters  after the descretization. For
instance the excitation strength `:A`:

```@example tutorial_09_forcing
rhs.params[:A]=1;
nothing #hide
```

!!! note
    You may also redefine the membrane impedance `:Y`, but note that you then
    also need to redefine it in `L` to get consistent results. There is
    currently no safety mechanism checking this for you!

## Solving the problem

As the problem is linear solving it is an easy matter. We just use the backslash
operator. For instance, finding the solution for an excitation frequency of 150 Hz
reads:

```@example tutorial_09_forcing
ω=150*2*pi;
sol=L(ω)\Array(rhs(ω));
nothing #hide
```

Note that the rhs is sparse, therefore, conversion to an Array is necessary as
the backslash operator does not support sparse rhs. Anyways, we have a solution!
We may right it to paraview, using the usual functions:

```@example tutorial_09_forcing
data=Dict();
data["150Hz abs"]=abs.(sol);
data["150Hz phase"]=angle.(sol);
vtk_write("forcing", mesh, data);
nothing #hide
```

## Summary

You can also discretize a source vector and use it to model forced problems.
Just use the `:speaker` option in the model description and set `soure=true`
in the call to discretize.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

