## standard example
using WavesAndEigenvalues.Helmholtz
mesh=Mesh("Rijke_mm.msh",scale=0.001)
dscrp=Dict() #initialize model discreptor
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
## refine mesh
fine_mesh=octosplit(mesh)
#solve fine mesh for comparison
c=generate_field(fine_mesh,speedofsound)
L=discretize(fine_mesh,dscrp,c)
sol,nn,flag=householder(L,400*2*pi,maxiter=10)
#write to file
data=Dict()
data["speed_of_sound"]=c
data["abs"]=abs.(sol.v)/maximum(abs.(sol.v)) #normalize so that max=1
data["phase"]=angle.(sol.v)
vtk_write("fine_test", fine_mesh, data) # Write the solution to paraview

## If you think that's not enough, well, then try...
fine_mesh=octosplit(fine_mesh)
