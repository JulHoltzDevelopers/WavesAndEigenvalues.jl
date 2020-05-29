"""
HD example and BC example
Dimensions:
tube length:    1000 mm
tube diameter:  20 mm
neck diameter:  5 mm
neck length:    20 mm
cavity diam:    20 mm
cavity len:     300 mm
Tube res freq = c/(2L) = 173.5 Hz
HD resonance = c/(2*pi)*sqrt(a_n/l_n*V_c) = 178.24 Hz
"""
##
using WavesAndEigenvalues.Helmholtz
using DelimitedFiles
include("fit_ss.jl")
##
mesh=Mesh("./WavesAndEigenvalues/examples/HD_full_geo.msh",scale=0.001) # Rescales the coordinates (all dymensions, with a mult)
C=ones(length(mesh.tetrahedra))*347.0
dscrp=Dict()
dscrp["Interior"]=(:interior,())

L=discretize(mesh, dscrp, C)
##
HD_data = readdlm("./WavesAndEigenvalues/examples/FR_HD.dat", '\t')
freq_HD = HD_data[:,1]
Z_HD = HD_data[:,2]+im*HD_data[:,3]
Npoles = 2
Ass, Bss, Css, Dss, poles, rmserr, fit = fit_ss_wrapper([1]./Z_HD,freq_HD,Npoles)

mesh1=Mesh("./WavesAndEigenvalues/examples/HD_eff_imp.msh",scale=0.001) # Rescales the coordinates (all dymensions, with a mult)
C1=ones(length(mesh1.tetrahedra))*347.0
dscrp1=Dict()
dscrp1["Interior"]=(:interior,())
dscrp1["ImpedanceBC"]=(:admittance,(Ass,Bss,Css,Dss))
#dscrp1["ImpedanceBC"]=(:admittance,(:Y,1e10))
L1=discretize(mesh1, dscrp1, C1)


##
sol, n, flag = householder(L,147*2*pi,maxiter=30, output = true, n_eig_val=1,tol=1e-9)
sol1, n1, flag1 = householder(L1,147*2*pi,maxiter=30, output = true, n_eig_val=1,tol=1e-9)

##
eigsBeyn, vecBeyn=beyn(L,[1 + 1im 2000+1im 2000-1im 1-1im], l=5, K=1, N=128,output=true,tol=1e-10)
eigsBeyn1, vecBeyn1=beyn(L1,[1 + 1im 2000+1im 2000-1im 1-1im], l=5, K=1, N=128,output=true)
