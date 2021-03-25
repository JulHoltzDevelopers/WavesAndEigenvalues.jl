```@meta
EditURL = "<unknown>/tutorial_06_second_order_elements.jl"
```

# Tutorial 06 Second Order Elements

So far we have deployed linear elements for the FEM discretization of the
Helmholtz equation. But WavesAndEigenvalues.Helmholtz can do more. Second order
elements are also available. All you have to do is to set the `order` keyword in
the call to `order=:quad`.

## Model set-up

We demonstrate things with the usual Rijke tube set-up from Tutorial 01.

```@example tutorial_06_second_order_elements
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
n=0.01 #interaction index
τ=0.001 #time delay
dscrp["Flame"]=(:flame,(γ,ρ,Q02U0,x_ref,n_ref,:n,:τ,n,τ)) #flame dynamics
R=287.05 # J/(kg*K) specific gas constant (air)
speedofsound(x,y,z) = z<0. ? sqrt(γ*R*Tu) : sqrt(γ*R*Tb)
c=generate_field(mesh,speedofsound)
```

The key difference is that we opt for second order elements during discretization

```@example tutorial_06_second_order_elements
L=discretize(mesh,dscrp,c, order=:quad)
```

All subsequent steps are the same. For instance, solving for an eigenvalue is done by

```@example tutorial_06_second_order_elements
sol,nn,flag=mslp(L,340*2*pi,maxiter=20,tol=1E-11,output=true)
#
```

And writing output to paraview by

```@example tutorial_06_second_order_elements
data=Dict()
data["abs mode"]=abs.(sol.v)
vtk_write("tutorial06",mesh, data)
```

Note However, that the vtk-file name will be `tutorial06_quad.vtu`.

## Summary
Quadratic elements are an easy matter. Just set `order=:quad` when calling
`discretize`.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

