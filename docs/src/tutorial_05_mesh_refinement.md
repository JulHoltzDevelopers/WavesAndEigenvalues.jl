```@meta
EditURL = "<unknown>/tutorial_05_mesh_refinement.jl"
```

# Tutorial 05 Mesh Refinement

## Intro
This tutorial explains you to utilities for performing mesh-refinement.
The main function here is `finer_mesh=octosplit(mesh)`, which splits all
tetrahedral cells in `mesh` by bisecting each edge. This will lead to a
subdivision of each tetrahedron into eight new tetrahedra. Such a division is
not unique. There are three ways of splitting a tetrahedron into eight new
tetrahedra and bisecting all six original edges. `octosplit` chooses the
splitting such that the maximum of the newly created edges is minimized.

## Basic set-up

We start with initializing the standard Rijke-Tube example.

```julia
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

```
Launching Householder...
Iter    Res:     dz:     z:
----------------------------------
0		Inf	Inf	2513.2741228718346
1		3.3678010621821997e6	672.6717790631509	1840.644381094955 + 7.520161245130288im
2		457309.5623269443	125.32039552244363	1715.3398137797765 + 9.511880065048567im
3		15754.920293139632	4.636826536115793	1710.7041295613465 + 9.614798251130926im
4		21.537680317163897	0.006356125798411251	1710.6977772512585 + 9.615018459466514im
5		4.045803703302238e-5	1.1939919720254716e-8	1710.6977772393393 + 9.615018460170372im
6		2.6334860638026723e-8	7.85077700900375e-12	1710.6977772393466 + 9.61501846017332im
7		1.4995064722670382e-8	4.324048085871753e-12	1710.6977772393423 + 9.615018460173136im
8		1.0838615959911618e-8	3.1913416131870233e-12	1710.697777239339 + 9.615018460173363im
9		1.4098455247391106e-8	4.099642246569255e-12	1710.6977772393432 + 9.615018460173125im
10		1.919724130587623e-9	4.655004020764636e-13	1710.6977772393427 + 9.615018460173225im
Warning: Maximum number of iterations has been reached!
...finished Householder!
#####################
 Householder results 
#####################
Number of steps: 10
Last step parameter variation:4.655004020764636e-13
Auxiliary eigenvalue λ residual (rhs):1.919724130587623e-9
Eigenvalue:1710.6977772393427 + 9.615018460173225im
Eigenvalue/(2*pi):272.26600738395945 + 1.5302777158563927im

```

## refine mesh
Now let's call `octosplit` for creating the refined mesh

```julia
fine_mesh=octosplit(mesh)
```

```
mesh: Rijke_mm.msh
#################
points:     6172
lines:      0
triangles:  6248
tetrahedra: 27040
#################
domains: Cold, Flame, Flame_in, Flame_out, Hot, Inlet, Interior, Outlet, Walls
```

and solve the problem again with the fine mesh for comparison.

```julia
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

```
mesh: Rijke_mm.msh
#################
points:     42507
lines:      0
triangles:  24992
tetrahedra: 216320
#################
domains: Cold, Flame, Flame_in, Flame_out, Hot, Inlet, Interior, Outlet, Walls
```

## Summary

The `octosplit`-function allows you to refine your mesh by bisecting all edges.
This creates a new mesh with eight-times as many tetrahedral cells. The tool is
quite useful for running h-convergence tests for your set-up.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

