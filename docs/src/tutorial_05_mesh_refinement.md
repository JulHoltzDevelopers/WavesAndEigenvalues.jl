```@meta
EditURL = "<unknown>/tutorial_05_mesh_refinement.jl"
```

# Tutorial 05 Mesh Refinement

## Intro
This tutorial introduces you to utilities for performing mesh-refinement.
The main comand here is `finer_mesh=octosplit(mesh)`, which splits all
tetrahedral cells in `mesh` by bisecting each edge. This will lead to a
subdivision of each tetrahedron into eight new tetrahedra. Such a division is
not unique. There are three ways of splitting a tetrahedron into eight new
tetrahedra and bisecting all six original edges. `octosplit` chooses the
splitting such that the maximum of the newly created edges is minimized.

## Basic set-up

We start with initializing the standard Rijke-Tube example.

```@example tutorial_05_mesh_refinement
using WavesAndEigenvalues.Helmholtz
mesh=Mesh("Rijke_mm.msh",scale=0.001)
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
sol,nn,flag=householder(L,400*2*pi,maxiter=10)
data=Dict()
data["speed_of_sound"]=c
data["abs"]=abs.(sol.v)/maximum(abs.(sol.v)) #normalize so that max=1
data["phase"]=angle.(sol.v)
vtk_write("test", mesh, data) # Write the solution to paraview
```

## refine mesh
Now let's call `octosplit` for creating the refined mesh

```@example tutorial_05_mesh_refinement
fine_mesh=octosplit(mesh)
```

and solve the problem again with the fine mesh for comparison.

```@example tutorial_05_mesh_refinement
c=generate_field(fine_mesh,speedofsound)
L=discretize(fine_mesh,dscrp,c)
sol,nn,flag=householder(L,400*2*pi,maxiter=10)
#write to file
data=Dict()
data["speed_of_sound"]=c
data["abs"]=abs.(sol.v)/maximum(abs.(sol.v)) #normalize so that max=1
data["phase"]=angle.(sol.v)
vtk_write("fine_test", fine_mesh, data) # Write the solution to paraview

# If you think that's not enough, well, then try...
fine_mesh=octosplit(fine_mesh)
```

## Summary

The `octosplit`-function allows you to refine your mesh by bisecting all edges.
This creates a new mesh with eight-times as many tetrahedral cells. The tool is
quite useful for running h-convergence tests for your set-up.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

