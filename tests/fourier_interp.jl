import WavesAndEigenvalues.NLEVP: generate_Σy_exp_ikx
using FFTW
##
N=7
B=ones(ComplexF64,N)
B[1]=0
Y=fft(B)/N
f=generate_Σy_exp_ikx(Y)
f(0.0+0.0im,0)
for z=0:N
    println("$z --> $(f(0.0+0.0im+z,0))")
end
