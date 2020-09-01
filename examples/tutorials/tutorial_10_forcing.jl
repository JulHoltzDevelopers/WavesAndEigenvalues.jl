# ## Model set-up.
# The model is the same Rijke tube configuration as in Tutorial 01:
using WavesAndEigenvalues.Helmholtz
mesh=Mesh("Rijke_mm.msh",scale=0.001) #load mesh
dscrp=Dict() #initialize model descriptor
dscrp["Interior"]=(:interior, ()) #define resonant cavity
dscrp["Outlet"]=(:speaker, (:A,1,:Y,1E15)) #specify outlet BC
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
L,rhs=discretize(mesh,dscrp,c,source=true)

rhs.params[:A]=1#E15*694/1im#694^2*1E15
## output
F=50:5:2000
tf=zeros(ComplexF64,size(F))
data=Dict()
for (idx,f)=enumerate(F)
    ω=2*pi*f
    mode="$(lpad(f,4,"0")) Hz"
    data[mode]=real.(L(ω)\Array(rhs(ω)))
    tf[idx]=data[mode][500]
end
#vtk_write("tut10", mesh, data)
##
#the speaker forcing is derived from an impedance boundary condition
#p̂ - (cZ)/(iω) * ∇p̂⋅n=A where p̂ is the pressure fluctuation amplitude, c the
#local speed of sound, Z a normalized impedance, i the imaginary unit, ω an
#angular frequency, n the outward pointing unit normal, and A the excitation
#level. This equation is equal to -iωY/c * p̂+∇p̂⋅n=-iωY/c*A, where the admittance
#is Y=1/Z. Note that for the special case of a a sound-soft impedance (Z=0 aka
#Y=∞) the excitation level A is identical to the pressure amplitude p̂ (p̂=A).
using PyPlot
###
plot(F,tf)

##
beyn
