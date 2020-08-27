
## https://julholtzdevelopers.github.io/WavesAndEigenvalues.jl/dev/tutorial_02_gloabal_eigenvalue_solver.html
using WavesAndEigenvalues.Helmholtz
mesh=Mesh("war_ne_geburt.msh",scale=0.001*0.75) #load mesh, scaling factor mm→m, scaling for 3D print
dscrp=Dict()
## create impedance model for mouth

H=0.015
W=0.004
R=0.0079
ΔL=2.3*R^2/sqrt(W*H)
ρ=1.225
Sp=pi*R^2
coeff=1im*ρ*ΔL/Sp
coeff/=ρ*347
coeff*=W*H
#coeff=1E20im
##
function admittance(z::ComplexF64,k::Int64)
    #f=coeff*z^-1
    f=factorial(k)
    if isodd(k)
        f*=-1
    end
    f*=z^(-1-k)
    f*=1/coeff #coefficient
    return f
end


##
dscrp["Interior"]=(:interior, ()) #define resonant cavity
# boundary conditions
dscrp["outlet"]=(:admittance, (:Y,1E15)) #specify outlet BC
#dscrp["mouth"]=(:admittance, (admittance,)) #specify mouth BC
dscrp["mouth"]=(:admittance, (:Y_mouth,.235E30)) #specify mouth BC
dscrp["inlet"]=(:admittance, (:Y,1E15)) #specify inlet BC
# the rest is rigid wall
speedofsound(x,y,z) = 347
c=generate_field(mesh,speedofsound) # creates array for each point in mesh according to speedofsound(x,y,z)
L=discretize(mesh,dscrp,c,el_type=1)
##
data=Dict()
sol,nn,flag=householder(L,1668*2*pi,maxiter=10,output=true)
data["speed_of_sound"]=c
data["abs:$(real(round(sol.params[:ω]/(2*pi)))) Hz"]=abs.(sol.v)/maximum(abs.(sol.v)) #normalize so that max=1
data["phase:$(real(round(sol.params[:ω]/(2*pi)))) Hz"]=angle.(sol.v)
##
#perturb_fast!(sol,L,:Y_mouth,30)
#sol(:Y_mouth,1000,15,15)/2/pi
##
# contour for eigenwert calculation
Γ=[20.0+5.0im, 20.0-5.0im, 15000.0-5.0im, 15000.0+5.0im].*2*pi
# wolf-juergen beyn eigenwert solver, Gauss-Legendre quadrature
Ω,P=beyn(L,Γ,N=32,l=100,output=true)
data=Dict()

##
for ω in Ω
sol,nn,flag=householder(L,ω,maxiter=10,output=false)
data["speed_of_sound"]=c
data["abs:$(real(round(sol.params[:ω]/(2*pi)))) Hz"]=abs.(sol.v)/maximum(abs.(sol.v)) #normalize so that max=1
data["phase:$(real(round(sol.params[:ω]/(2*pi)))) Hz"]=angle.(sol.v)
end

# data["speed_of_sound"]=c
# data["abs:$(real(round(sol.params[:ω]/(2*pi)))) Hz"]=abs.(sol.v)/maximum(abs.(sol.v)) #normalize so that max=1
# data["phase:$(real(round(sol.params[:ω]/(2*pi)))) Hz"]=angle.(sol.v)
vtk_write("war_ne_geburt", mesh, data) # Write the solution to paraview
