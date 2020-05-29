import LinearAlgebra: norm, qr, diagm, eigvals

function fit_ss_wrapper(f,freqs,Npoles)
    s = freqs*2*pi*im #Convert frequency in Lalplace variable

    poles = Array{ComplexF64}(undef,0)
    if mod(Npoles,2)==1
        push!(poles,-1)
        for ii=1:(Npoles-1)/2
            push!(poles, -1-im*2*pi*maximum(freqs)/ii, -1+im*2*pi*maximum(freqs)/ii)
        end
    else
        for ii=1:Npoles/2
            push!(poles, -1-im*2*pi*maximum(freqs)/ii, -1+im*2*pi*maximum(freqs)/ii)
        end
    end

A, B, C, D, poles, rmserr, fit = fit_ss(f,s,poles)

# Call multiple times to ensure convergence of poles
for i=1:9
    A, B, C, D, poles, rmserr, fit = fit_ss(f,s,poles)
end

return A, B, C, D, poles, rmserr, fit

end

# Code adapted from vectorfit
function fit_ss(f,s,poles)
# f: frequency response to fit
# s: Laplace variable=sigma+im*omega, at which f has been evaluated
# poles: guess of the location of the poles.
#Note: The length(poles) is the order of the ss approximation

#Returns
#A, B, C, D: state-space matrices)
#fit: the frequency response is reconctructed as C*(s*I-A)^{-1}*B + D
#rmserr: rms error between data nad fit
#poles: the location of the converged poles (Note: run the routine multiple times to ensure convergence)

if size(s,1)<size(s,2)
    s=transpose(s)
end

if size(poles,2)<size(poles,1)
    poles=transpose(poles)
end

if size(f,2)<size(f,1)
    f=transpose(f)
end

# Should include here some sanity checks
TOLlow=1e-18
TOLhigh=1e18
LAMBD=diagm(transpose(poles))
Ns=length(s)
N=size(LAMBD,1)
Nc=length(f[:,1])
weight = transpose(ones(size(f)))
B=ones(N,1)
SERA=poles
SERC=zeros(Nc,N)
SERD=zeros(Nc,1)
offs = 0
roetter=poles
fit=zeros(ComplexF64,Nc,Ns)

#################################################
# Identify poles

Escale=zeros(1,Nc+1)

# Finding out which starting poles are complex:
cindex=zeros(1,N)
for m=1:N
    if imag(LAMBD[m,m])!=0
        if m==1
            cindex[m]=1
        else
            if cindex[m-1]==0 || cindex[m-1]==2
                cindex[m]=1
                cindex[m+1]=2
            else
                cindex[m]=2
            end
        end
    end
end

#################################################
# Building system - matrix :

Dk=zeros(ComplexF64,Ns,N+1)
for m=1:N
    if cindex[m]==0      #real pole
        Dk[:,m]=[1]./(s.-LAMBD[m,m])
    elseif cindex[m]==1  #complex pole, 1st part
        Dk[:,m]  =[1]./(s.-LAMBD[m,m]) + [1]./(s.-LAMBD[m,m]')
        Dk[:,m+1]=[im]./(s.-LAMBD[m,m]) - [im]./(s.-LAMBD[m,m]')
    end
end

Dk[:,N+1].=1

#Scaling for last row of LS-problem (pole identification)
scale=0
for m=1:Nc
    if length(weight[1,:])==1
        scale=scale+(norm(weight[:,1].*f[m,:]))^2
    else
        scale=scale+(norm(weight[:,m].*f[m,:]))^2
    end
end
scale=sqrt(scale)/Ns

AA=zeros(Nc*(N+1),N+1)
bb=zeros(Nc*(N+1),1)
Escale=zeros(1,length(AA[1,:]))
for n=1:Nc
    A=zeros(ComplexF64,Ns,(N+offs) +N+1)
    weig=weight[:,n]

    for m=1:N+offs #left block
        A[1:Ns,m]=weig.*Dk[1:Ns,m]
    end
    inda=N+offs
    for m=1:N+1 #right block
        A[1:Ns,inda+m]=-weig.*Dk[1:Ns,m].*f[n,1:Ns]
    end

    A=[real(A);imag(A)]

    #Integral criterion for sigma:
    offset=(N+offs)
    temp=zeros(1,size(A,2))
    if n==Nc
        for mm=1:N+1
            temp[offset+mm]=real(scale*sum(Dk[:,mm]))
        end
    end
    A = [A;temp]

    QR_A=qr(A) #AO: this would be more efficient with the reduced QR decomposition
    Q = QR_A.Q
    R = QR_A.R
    ind1=N+offs+1
    ind2=N+offs+N+1
    R22=R[ind1:ind2,ind1:ind2]
    AA[(n-1)*(N+1)+1:n*(N+1),:]=R22
    if n==Nc
        bb[(n-1)*(N+1)+1:n*(N+1),1]=Q[end,N+offs+1:N+offs+1+N]'*Ns*scale
    end
end

for col=1:length(AA[1,:])
    Escale[col]=1/norm(AA[:,col])
    AA[:,col]=Escale[col].*AA[:,col]
end

x=AA\bb
x=x.*transpose(Escale)

#if produced D of sigma extremely small and large. Solve again, without relaxation
if abs(x[end])<TOLlow || abs(x[end])>TOLhigh
    AA=zeros(Nc*(N),N)
    bb=zeros(Nc*(N),1)

    if x[end]==0
        Dnew=1
    elseif abs(x[end])<TOLlow
        Dnew=sign(x[end])*TOLlow
    elseif abs(x[end])>TOLhigh
        Dnew=sign(x[end])*TOLhigh
    end

    Escale=zeros(1,N)

    for n=1:Nc
        A=zeros(ComplexF64,Ns,(N+offs) +N)
        weig=weight[:,n]

        for m=1:N+offs #left block
            A[1:Ns,m]=weig.*Dk[1:Ns,m]
        end
        inda=N+offs
        for m=1:N #right block
            A[1:Ns,inda+m]=-weig.*Dk[1:Ns,m].*f[n,1:Ns]
        end

        b=Dnew*weig.*f[n,1:Ns]
        A=[real(A);imag(A)]
        b=[real(b);imag(b)]
        offset=(N+offs)
        QR_A=qr(A) #AO: this would be more efficient with the reduced QR decomposition
        Q = QR_A.Q
        R = QR_A.R
        ind1=N+offs+1
        ind2=N+offs+N
        R22=R[ind1:ind2,ind1:ind2]
        AA[(n-1)*N+1:n*N,:]=R22
        bb[(n-1)*N+1:n*N,1]=transpose(Q[:,ind1:ind2])*b
    end

    for col=1:length(AA[1,:])
        Escale[col]=1/norm(AA[:,col])
        AA[:,col]=Escale[col].*AA[:,col]
    end

    x=AA\bb
    x=x.*transpose(Escale)
    x=[x;Dnew]
end

C=x[1:end-1]
C = convert(Array{ComplexF64,1},C)
D=x[end]

#################################################

#We now change back to make C complex :
for m=1:N
    if cindex[m]==1
        for n=1:1%Nc+1
            r1=C[m]
            r2=C[m+1]
            C[m]=r1+im*r2
            C[m+1]=r1-im*r2
        end
    end
end

#################################################
# We now calculate the zeros for sigma: This converts ALL to be real!
m=0
for n=1:N
    m=m+1
    if m<N
        if( abs(LAMBD[m,m])>abs(real(LAMBD[m,m])) )
            LAMBD[m+1,m]=-imag(LAMBD[m,m])
            LAMBD[m,m+1]=imag(LAMBD[m,m])
            LAMBD[m,m]=real(LAMBD[m,m])
            LAMBD[m+1,m+1]=LAMBD[m,m]
            B[m,1]=2
            B[m+1,1]=0
            koko=C[m]
            C[m]=real(koko)
            C[m+1]=imag(koko)
            m=m+1
        end
    end
end

ZER=real(LAMBD-B*transpose(C)/D) #AO: convert to real! This ensures that eigenvalues are either real or complex-conjugate
roetter=transpose(eigvals(ZER))
unstables=real(roetter).>0
#Forcing unstable poles to be stable...
roetter[unstables]=roetter[unstables]-2*real(roetter[unstables])

roetter=transpose(sort(roetter[1,:],by=abs))
N=length(roetter)


#################################################
#sorte poles s.t. the real come first:
for n=1:N
    for m=n+1:N
        if imag(roetter[m])==0 && imag(roetter[n])!=0
            trans=roetter[n]
            roetter[n]=roetter[m]
            roetter[m]=trans
        end
    end
end

N1=0
for m=1:N
    if imag(roetter[m])==0
         N1=m
      end
end
if N1<N
    roetter[N1+1:N]=transpose(sort(roetter[N1+1:N],by=abs))
end

roetter=roetter-2*im*imag(roetter)
SERA=transpose(roetter)

#################################################
# RESIDUE IDENTIFICATION:
# We now calculate SER for f, using the modified zeros of sigma as new poles :

LAMBD=roetter
# Finding out which poles are complex :
cindex=zeros(1,N);
for m=1:N
    if imag(LAMBD[m])!=0
        if m==1
            cindex[m]=1
        else
            if cindex[m-1]==0 || cindex[m-1]==2
                cindex[m]=1
                cindex[m+1]=2
            else
                cindex[m]=2
            end
        end
    end
end


#################################################
# We now calculate the SER for f (new fitting), using the above calculated
# zeros as known poles :

A=zeros(ComplexF64,2*Ns,N)
BB=zeros(ComplexF64,2*Ns,Nc)
Dk=zeros(ComplexF64,Ns,N);

for m=1:N
  if cindex[m]==0      #real pole
    Dk[:,m]=weight./(s.-LAMBD[m])
  elseif cindex[m]==1  #complex pole, 1st part
    Dk[:,m]  =weight./(s.-LAMBD[m]) + weight./(s.-LAMBD[m]')
    Dk[:,m+1]=im.*weight./(s.-LAMBD[m]) - im.*weight./(s.-LAMBD[m]')
  end
end

A[1:Ns,1:N]=Dk

for m=1:Nc
    BB[1:Ns,m]=weight.*f[m,:]
end

A[Ns+1:2*Ns,:]=imag(A[1:Ns,:])
A[1:Ns,:]=real(A[1:Ns,:])
BB[Ns+1:2*Ns,:]=imag(BB[1:Ns,:])
BB[1:Ns,:]=real(BB[1:Ns,:])

Escale=zeros(1,length(A[1,:]))
for col=1:length(A[1,:])
    Escale[col]=norm(A[:,col],2)
    A[:,col]=A[:,col]./Escale[col]
end

X=A\BB;

for n=1:Nc
    X[:,n]=X[:,n]./transpose(Escale)
end

X=transpose(X)
C=X[:,1:N]


#################################################
#We now change back to make C complex.

for m=1:N
  if cindex[m]==1
    for n=1:Nc
      r1=C[n,m]
      r2=C[n,m+1]
      C[n,m]=r1+im*r2
      C[n,m+1]=r1-im*r2
   end
  end
end

B=ones(N,1);

SERA  = LAMBD;
SERB  = B;
SERC  = C;

Dk=zeros(ComplexF64,Ns,N)
for m=1:N
    Dk[:,m]=[1]./(s.-SERA[m])
end
for n=1:Nc
    fit[n,:]=transpose((Dk*SERC[n,:]))
end

fit=transpose(fit)
f=transpose(f)
diff=fit-f
rmserr=sqrt(sum(sum(abs.(diff.^2))))/sqrt(Nc*Ns)
fit=transpose(fit)

A=SERA
poles=A
B=SERB
C=SERC
D=SERD

#################################################
# Convert into real state-space model

A=diagm(transpose(A))
cindex=zeros(1,N)

for m=1:N
if imag(A[m,m])!=0
  if m==1
    cindex[m]=1
  else
    if cindex[m-1]==0 || cindex[m-1]==2
      cindex[m]=1
       cindex[m+1]=2
    else
      cindex[m]=2
    end
  end
end
end

n=0;
for m=1:N
    n=n+1;
    if cindex[m]==1
        a=A[n,n];
        a1=real(a);
        a2=imag(a);
        c=C[:,n];
        c1=real(c);
        c2=imag(c);
        b=B[n,:];
        b1=2*real(b);
        b2=-2*imag(b);
        Ablock=[a1 a2;-a2 a1];

        A[n:n+1,n:n+1]=Ablock;
        C[:,n]=c1;
        C[:,n+1]=c2;
        B[n,:]=b1;
        B[n+1,:]=b2;
    end
end


return A, B, C, D, poles, rmserr, fit


end
