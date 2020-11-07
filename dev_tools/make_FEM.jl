import SymPy
x, y, z,a = SymPy.symbols("x y z a")
b1, b2, b3, b4, b5, b6, b7, b8, b9, b10 = SymPy.symbols("b1 b2 b3 b4 b5 b6 b7 b8 b9 b10")
c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 = SymPy.symbols("c1 c2 c3 c4 c5 c6 c7 c8 c9 c10")
A11, A12, A13, A22, A23, A33 = SymPy.symbols("A11 A12 A13 A22 A23 A33")
A31, A21, A32 = A13, A12, A23
B11, B12, B13, B21, B22, B23, B31, B32, B33 = SymPy.symbols("B11 B12 B13 B21 B22 B23 B31 B32 B33")


A=SymPy.Matrix([A11 A12 A13;
              A21 A22 A23;
              A31 A32 A33])
B=SymPy.Matrix([B11 B12 B13;
                B21 B22 B23;
                B31 B32 B33])

a=1-x-y-z


#second order shape functions
f1 = (2*x-1)*x
f2 = (2*y-1)*y
f3 = (2*z-1)*z
f4 = (2*a-1)*a
f5 = 4*x*y
f6 = 4*x*z
f7 = 4*x*a
f8 = 4*y*z
f9 = 4*y*a
f10 = 4*z*a
##
#hermite shape functions (third order polynomial)
fh_vtx(x,y,z)=1-3*x^2-13*x*y-13*x*z-3*y^2-13*y*z-3*z^2+2*x^3+13*x^2*y+13*x^2*z+13*x*y^2+33*x*y*z+13*x*z^2+2*y^3+13*y^2*z+13*y*z^2+2*z^3
fh_dvtx(x,y,z)=x-2*x^2-3*x*y-3*x*z+x^3+3*x^2*y+3*x^2*z+2*x*y^2+4*x*y*z+2*x*z^2
fh_dvtxx(x,y,z)=-x^2+2*x*y+2x*z+x^3-2*x^2*y-2*x^2*z-2*x*y^2-2x*y*z-2*x*z^2
fh_fc(x,y,z)=27*x*y*z
##
fh01=fh_vtx(y,z,a)
fh02=fh_vtx(x,z,a)
fh03=fh_vtx(x,y,a)
fh04=fh_vtx(x,y,z)
#partial derivatives (arguement order matters)
# wrt x
fh05=fh_dvtxx(x,y,z)
fh06=fh_dvtx(x,a,z).expand()
fh07=fh_dvtx(x,y,a).expand()
fh08=fh_dvtx(x,y,z).expand()
# wrt y
fh09=fh_dvtx(y,z,a).expand()
fh10=fh_dvtxx(y,x,z)
fh11=fh_dvtx(y,x,a).expand()
fh12=fh_dvtx(y,x,z).expand()

#wrt z
fh13=fh_dvtx(z,y,a).expand()
fh14=fh_dvtx(z,x,a).expand()
fh15=fh_dvtxx(z,x,y)
fh16=fh_dvtx(z,x,y).expand()
#
fh17=fh_fc(y,z,a)
fh18=fh_fc(x,z,a)
fh19=fh_fc(x,y,a)
fh20=fh_fc(x,y,z)

#recombined implementation

fr05=(B11*fh05+B12*fh09+B13*fh13).simplify()
fr06=(B11*fh06+B12*fh10+B13*fh14).simplify()
fr07=(B11*fh07+B12*fh11+B13*fh15).simplify()
fr08=(B11*fh08+B12*fh12+B13*fh16).simplify()

fr09=(B21*fh05+B22*fh09+B23*fh13).simplify()
fr10=(B21*fh06+B22*fh10+B23*fh14).simplify()
fr11=(B21*fh07+B22*fh11+B23*fh15).simplify()
fr12=(B21*fh08+B22*fh12+B23*fh16).simplify()

fr13=(B31*fh05+B32*fh09+B33*fh13).simplify()
fr14=(B31*fh06+B32*fh10+B33*fh14).simplify()
fr15=(B31*fh07+B32*fh11+B33*fh15).simplify()
fr16=(B31*fh08+B32*fh12+B33*fh16).simplify()




##
function check_hermitian(f)
    for (j,pnt) in enumerate([(1,0,0), (0,1,0), (0,0,1), (0,0,0),(0,1/3,1/3), (1/3,0,1/3), (1/3,1/3,0), (1/3,1/3,1/3)])
        xx,yy,zz=pnt
        res=f.subs(x,xx).subs(y,yy).subs(z,zz).simplify()
        println("$pnt ---> $res")
        res=diff(f,x).subs(x,xx).subs(y,yy).subs(z,zz).simplify()
        if j<=4
            println("dx $pnt ---> $res")
            res=diff(f,y).subs(x,xx).subs(y,yy).subs(z,zz).simplify()
            println("dy $pnt ---> $res")
            res=diff(f,z).subs(x,xx).subs(y,yy).subs(z,zz).simplify()
            println("dz $pnt ---> $res")
        end
    end
end


##
#zweiter ordnung
X=[x,y,z,a]
F= [[x,y,z,a], [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10 ],
    [fh01, fh02,fh03,fh04,fh05,fh06,fh07,fh08,fh09,fh10,fh11,fh12,fh13,fh14,fh15,fh16,fh17,fh18,fh19,fh20],
    [fh01, fh02,fh03,fh04,fr05,fr06,fr07,fr08,fr09,fr10,fr11,fr12,fr13,fr14,fr15,fr16,fh17,fh18,fh19,fh20]]
C= [[c1*x, c2*y, c3*z, c4*a,],[c1*f1, c2*f2, c3*f3, c4*f4, c5*f5, c6*f6, c7*f7, c8*f8, c9*f9, c10*f10 ],]
b= [[b1*x, b2*y, b3*z, b4*a,],[b1*f1, b2*f2, b3*f3, b4*f4, b5*f5, b6*f6, b7*f7, b8*f8, b9*f9, b10*f10 ],]

f = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10 ]
ff=[f1,f4,f5,f7]
surf=[[1,2,4],[1,2,4,5,7,9],[1 2 4 5 6 8 9 10 12 13 14 16 19],[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]]


##
function check_hermitian_mat(F)
    M=Array{Any}(undef,20,20)
    for (i,f) in enumerate(F)
        for (j,pnt) in enumerate([(1,0,0), (0,1,0), (0,0,1), (0,0,0), (0,1/3,1/3), (1/3,0,1/3), (1/3,1/3,0), (1/3,1/3,1/3)])
            xx,yy,zz=pnt
            res=f.subs(x,xx).subs(y,yy).subs(z,zz).simplify()
            #println("$pnt ---> $res")

            res=diff(f,x).subs(x,xx).subs(y,yy).subs(z,zz).simplify()
            if j >4
                idx=j+12
            else
                idx=j
            end
            M[idx,i]=res

            if j<=4
                res=diff(f,x).subs(x,xx).subs(y,yy).subs(z,zz).simplify()
                #println("dx $pnt ---> $res")
                M[j+4,i]=res
                res=diff(f,y).subs(x,xx).subs(y,yy).subs(z,zz).simplify()
                #println("dy $pnt ---> $res")
                M[j+8,i]=res
                res=diff(f,z).subs(x,xx).subs(y,yy).subs(z,zz).simplify()
                #println("dz $pnt ---> $res")
                M[j+12,i]=res
            end
        end
    end
    return M
end
TT=check_hermitian_mat(F[3])
#erster ordnung
#f = [x y z a]
# cs=[c1*x, c2*y, c3*z, c4*a,]
# c=sum(cs)
# #0. ordnung
# c=SymPy.symbols("c")
##

function grad(f)
    return [SymPy.diff(f,x); SymPy.diff(f,y); SymPy.diff(f,z)]
end

function matmul(h,A,g)
    out=0
    for i = 1:3
        for j=1:3
            out+=h[i]*A[i,j]*g[j]
        end
    end
    return out
end
##
function gradient(f)
    N=length(f)
    return transpose(SymPy.Matrix([SymPy.diff(f,x) SymPy.diff(f,y) SymPy.diff(f,z)]))
end

function build_gradient_matrix(f)
    M=gradient(f[0])
    for i = 1:N-1
        M=M.col_insert(i,gradient(f[i]).T)
    end
    return M
end

function tet_integrate(f)
    if size(f)!=()
        f=f[1]
    end
    return SymPy.integrate(SymPy.integrate(SymPy.integrate(f, (x, 0, 1 - z -y)), (y, 0, 1 - z)), (z, 0, 1))
end

function tri_integrate(f)
    return SymPy.integrate(SymPy.integrate(f, (x, 0, 1 - y)), (y, 0, 1))
end


function mass_matrix(f_i,f_j)
    #f_i=f_i.subs(a,suba)
    #f_j=f_j.subs(a,suba)
    return tet_integrate(f_i*f_j)
end

function stiff_matrix(f_i,f_j)
    f_i=gradient(f_i)
    f_j=gradient(f_j)
    #f_i=f_i.subs(a,suba)
    #f_j=f_j.subs(a,suba)
    return -tet_integrate(c^2*transpose(f_i)*A*f_j)
end

function source_vector(f)
    #f=f.subs(a,suba)
    return tet_integrate(f)
end

function boundary_mass_matrix(f_i, f_j)
    f = f_i*f_j
    return tri_integrate(c*f).subs(z, 0)
end

function boundary_source_vector(f)
    #f=f.subs(a,suba)
    return tri_integrate(c*f).subs(z, 0)
end

## build (q,p_j)
function local_matrix(order=1,corder=2, crdnt=0, crdnt_adj=0 ;symmetric=false,typ=:none)
    f=F[order]
    if typeof(corder)==SymPy.Sym
        coeff=corder
    elseif corder==0
        coeff=1
    else
        coeff=sum(C[corder])
    end

    if typ==:boundary
        f=f[surf[order]]
    end

    M=Array{Any}(undef,length(f),length(f))
        for (i,fi) in enumerate(f)
            if typ==:nabla
                    fi=grad(fi)
            end
            for (j,fj) in enumerate(f)
                if symmetric && j<i
                    continue
                end
                println("$i $j"); flush(stdout)
                if typ==:nabla
                    fj=grad(fj)
                    intg=matmul(fj,A,fi)
                elseif typ==:boundary
                    intg=(fj.subs(z,0)*fi.subs(z,0)).simplify()
                elseif length(crdnt) == 2
                    intg=fi*SymPy.diff(fj, X[crdnt[1]], X[crdnt[2]])
                elseif  crdnt==0 && crdnt_adj==0
                    intg=(fi*fj).simplify()
                elseif crdnt!=0 && crdnt_adj==0
                    intg=(fi*SymPy.diff(fj,X[crdnt])).simplify()
                elseif crdnt==0 && crdnt_adj!=0
                    intg=(SymPy.diff(fi,X[crdnt_adj])*fj).simplify()
                elseif crdnt!=0 && crdnt_adj!=0
                    intg=(SymPy.diff(fi,X[crdnt_adj])*SymPy.diff(fj,X[crdnt])).simplify()
                end
                if typ==:boundary
                    M[i,j]=tri_integrate((coeff*intg).subs(z,0))
                else
                    M[i,j]=tet_integrate(coeff*intg)
                end
                if symmetric
                    M[j,i]=M[i,j]
                end
            end
        end

    return M
end

function local_vector(order=1,corder=2, crdnt=0, crdnt_adj=0, typ=:none )
    f=F[order]
    if typeof(corder)==SymPy.Sym
        coeff=corder
    elseif corder==0
        coeff=1
    else
        coeff=sum(C[corder])
    end
    if typ==:boundary
        f=f[surf[order]]
    end

    M=Array{Any}(undef,length(f))
    for (i,fi) in enumerate(f)
        if typ==:boundary
            M[i]=tri_integrate((fi*coeff).subs(z, 0).simplify()).simplify()
        else
            M[i]=tet_integrate(fi*coeff).simplify()
        end
    end

    return M
end

##

#Dx=local_matrix(1,diff(sum(C[1])),1)
#Dy=local_matrix(1,0,2)
#Dz=local_matrix(1,0,3)

##
N=4
function print_matrix(M,coeff=false,symmetric=false;diff="")
    N=size(M,1)
    if coeff
        txt=""
        for i=1:N
            for j=(symmetric ? i : 1) : N
                txt*="M$diff[$i,$j]="*string(M[i,j])*"\n"
            end
        end
        if symmetric
            for i=1:N
                for j=i+1:N
                    txt*="M$diff[$j,$i]=M$diff[$i,$j]\n"
                end
            end
        end
    else
        txt="["
        for i=1:N
            for j=1:N
                if j!=N
                    txt*=string(M[i,j])*" "
                else
                    txt*=string(M[i,j])
                end
            end
            if i!=N
                txt*=";\n"
            else
                txt*="]"
            end
        end
    end

    txt=replace(txt,"A11"=>"A[1,1]")
    txt=replace(txt,"A12"=>"A[1,2]")
    txt=replace(txt,"A13"=>"A[1,3]")
    txt=replace(txt,"A21"=>"A[2,1]")
    txt=replace(txt,"A22"=>"A[2,2]")
    txt=replace(txt,"A23"=>"A[2,3]")
    txt=replace(txt,"A31"=>"A[3,1]")
    txt=replace(txt,"A32"=>"A[3,2]")
    txt=replace(txt,"A33"=>"A[3,3]")

    txt=replace(txt,"B11"=>"B[1,1]")
    txt=replace(txt,"B12"=>"B[1,2]")
    txt=replace(txt,"B13"=>"B[1,3]")
    txt=replace(txt,"B21"=>"B[2,1]")
    txt=replace(txt,"B22"=>"B[2,2]")
    txt=replace(txt,"B23"=>"B[2,3]")
    txt=replace(txt,"B31"=>"B[3,1]")
    txt=replace(txt,"B32"=>"B[3,2]")
    txt=replace(txt,"B33"=>"B[3,3]")

    for i=1:10
        txt=replace(txt,"c$i^2"=>"cc[$i,$i]")
        txt=replace(txt,"b$i^2"=>"bb[$i,$i]")
        for j=1:10
            txt=replace(txt,"c$i*c$j"=>"cc[$i,$j]")
            txt=replace(txt,"b$i*c$j"=>"bc[$i,$j]")
            txt=replace(txt,"b$i*b$j"=>"bb[$i,$j]")
        end
    end


    println(txt)
    return txt
end
##
#M =local_matrix(1,0) mass matrix
#M=local_matrix(1,0,0,3)
#M=local_matrix(1,diff(sum(C[1]),z),0)
for i=1:3
    M=local_matrix(3,sum(b[2])*sum(C[1]),i,typ=:diff,symmetric=false)
    txt=print_matrix(M,true,false,diff="$i")
    open("M$i$i.jl","w") do fid
        write(fid,txt)
    end
end
## vectors

M=local_vector(3,0,0,0,:boundary)
#M[herm_surf]
###
### derivatives
function deriv()
    DD=Array{Any}(undef,4,3)
    df=0
    for j=1:3
        for fi in C[2]
            for (i,xx) in enumerate([x,y,z])
                df+=SymPy.diff(fi,x)*A[i,j]
            end
        end
        DD[1,j]=dfi.subs(x,1).subs(y,0).subs(z,0)
        DD[2,j]=dfi.subs(x,0).subs(y,1).subs(z,0)
        DD[3,j]=dfi.subs(x,0).subs(y,0).subs(z,1)
        DD[4,j]=dfi.subs(x,0).subs(y,0).subs(z,0)
    end
    return DD
end

##
DD=deriv()
###############################################################
## legacy code
###############################################################
N=length(f)

M=Array{Any}(undef,N,N)
K=Array{Any}(undef,N,N)
C=Array{Any}(undef,N,N)
S=Array{Any}(undef,N)
G=Array{Any}(undef,N)
H=Array{Any}(undef,N)

for (i,fi) in enumerate(f)
    S[i] = source_vector(fi)
    G[i] = gradient(fi)
    H[i] =  boundary_source_vector(fi)
    for (j,fj) in enumerate(f)
        if j<i
            continue
        end
        println("$i $j"); flush(stdout)
        M[i,j]=mass_matrix(fi,fj)
        K[i,j]=stiff_matrix(fi,fj)
        C[i,j]=boundary_mass_matrix(fi,fj)
    end
end

##
function print_K(K)
    I,J=size(K)
    out=""
    for i =1:I
        for j= 1:J
            if i<=j
                txt=string(K[i,j])
                txt=replace(txt,"A11"=>"A[1,1]")
                txt=replace(txt,"A12"=>"A[1,2]")
                txt=replace(txt,"A13"=>"A[1,3]")
                txt=replace(txt,"A21"=>"A[2,1]")
                txt=replace(txt,"A22"=>"A[2,2]")
                txt=replace(txt,"A23"=>"A[2,3]")
                txt=replace(txt,"A31"=>"A[3,1]")
                txt=replace(txt,"A32"=>"A[3,2]")
                txt=replace(txt,"A33"=>"A[3,3]")
            else
                txt="m[$j,$i]"
            end
            out*="m[$i,$j]="*txt*"\n"
            #println("m[$i,$j]="*txt)
        end
    end

    return out
end
##
txt=print_K(K)
## quad prints
N=10
function print_M(M)
    out="["
    txt=""
    for i=1:N
        for j=1:N
            if i<=j
                txt*=string(M[i,j])
            else
                txt*=string(M[j,i])
            end
            txt*=" "
        end
        txt*=";\n"
    end
    return txt
end

txtM=print_M(M)

##
function print_K(K)
    txt=""
    for i=1:N
        for j=1:N
            if i<=j
                temp=string(K[i,j])
                temp=replace(temp,"A11"=>"A[1,1]")
                temp=replace(temp,"A12"=>"A[1,2]")
                temp=replace(temp,"A13"=>"A[1,3]")
                temp=replace(temp,"A21"=>"A[2,1]")
                temp=replace(temp,"A22"=>"A[2,2]")
                temp=replace(temp,"A23"=>"A[2,3]")
                temp=replace(temp,"A31"=>"A[3,1]")
                temp=replace(temp,"A32"=>"A[3,2]")
                temp=replace(temp,"A33"=>"A[3,3]")
                temp=replace(temp,"c1^2"=>"c[1,1]")
                temp=replace(temp,"c2^2"=>"c[2,2]")
                temp=replace(temp,"c3^2"=>"c[3,3]")
                temp=replace(temp,"c4^2"=>"c[4,4]")
                temp=replace(temp,"c1*c2"=>"c[1,2]")
                temp=replace(temp,"c1*c3"=>"c[1,3]")
                temp=replace(temp,"c1*c4"=>"c[1,4]")
                temp=replace(temp,"c2*c3"=>"c[2,3]")
                temp=replace(temp,"c2*c4"=>"c[2,4]")
                temp=replace(temp,"c3*c4"=>"c[3,4]")

                txt*="m[$i,$j]="*temp
                txt*="\n"
            else
                #txt*="m[$i,$j]=m[$j,$i]"
            end

        end
    end
    return txt
end





txtK=print_K(K)
##
function print_C(C)
    D=Dict()
    D[1]=1
    D[2]=2
    D[4]=3
    D[5]=4
    D[7]=5
    D[9]=6
    txt=""
    for i=1:N
        for j=1:N
            if i in [3,6,8,10] || j in [3,6,8,10]
                continue
            end
            I=D[i]
            J=D[j]

            if i<=j
                temp=string(C[i,j])
                temp=replace(temp,"c4"=>"c3")
                txt*="m[$I,$J]="*temp
                txt*="\n"
            else
                #txt*="m[$i,$j]=m[$j,$i]"
            end

        end
    end
    return txt
end


txtC=print_C(C)
println(txtC)
##
function print_G(G)
    txt=""
    for i=1:10
        for j =1:3
            temp=string(G[i][j])
            temp=replace(temp," "=>"")
        txt*=temp*" "
        end
        txt*=";\n"
    end
    return txt
end
 txtG=print_G(G)
 println(txtG)

##
 for i=1:10
     for j=1:10
         if j>=i
             continue
         end
        println("m[$i,$j]=m[$j,$i]")
    end
 end

##shape sensitivity
 #first term
q1, q2, q3, q4, q5, q6, q7, q8, q9, q10=SymPy.symbols("q1 q2 q3 q4 q5 q6 q7 q8 q9 q10")
p1, p2, p3, p4, p5, p6, p7, p8, p9, p10=SymPy.symbols("p1 p2 p3 p4 p5 p6 p7 p8 p9 p10")
q=[q1, q2, q3, q4, q5, q6, q7, q8, q9, q10]
p=[p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]
n1, n2, n3=SymPy.symbols("n1 n2 n3")
#p_adj=q1*f1+q2*f2+q4*f4+q5*f5+q7*f7+q9*f9
#p_adj=q1*f1+q2*f2+q3*f3+q4*f4+q5*f5+q6*f6+q7*f7+q8*f8+q9*f9+q10*f10
#p=p1*f1+p2*f2+p3*f3+p4*f4+p5*f5+p6*f6+p7*f7+p8*f8+p9*f9+p10*f10
cc=SymPy.symbols("c")
ξ=[x,y,z,a]
##
i=j=k=l=1
#TODO: do the summations
##Term I
result=tri_integrate.(f[j]*gradient(f[i])*ξ[k])
##Term II
res=gradient(f[i]).*ξ[k]
term=[]
for func in res
    append!(term,f[j].*gradient(func))
end
erm=reshape(term,3,3)
term=tri_integrate.(term)
##TermIII
term3=tri_integrate.(f[j].*gradient(f[i]).*ξ[k])



#############################################
#Legacy Code

##res=tri_integrate(cc^2*p_adj*transpose(gradient(p)))
res=tri_integrate(cc^2*q1*f1*transpose(gradient(p1*f1)))*0
for (j,(qq,fq)) in enumerate(zip(q,f))
    for (i,(pp,fp)) in enumerate(zip(p,f))
        println("$j:$i")
        global res
        tmp=fq.subs(z,-z)*transpose(gradient(fp.subs(z,-z)))*a
        for idx=1:3
            tmp[idx]=qq*tmp[idx].subs(z,0)*pp
        end
        res+=tri_integrate(tmp)
    end
end
##
function print_term(term)
    temp=string(term)
    temp=replace(temp,"c1^2"=>"cc[1,1]")
    temp=replace(temp,"c2^2"=>"cc[2,2]")
    temp=replace(temp,"c3^2"=>"cc[3,3]")
    temp=replace(temp,"c4^2"=>"cc[4,4]")
    temp=replace(temp,"c1*c2"=>"cc[1,2]")
    temp=replace(temp,"c1*c3"=>"cc[1,3]")
    temp=replace(temp,"c1*c4"=>"cc[1,4]")
    temp=replace(temp,"c2*c3"=>"cc[2,3]")
    temp=replace(temp,"c2*c4"=>"cc[2,4]")
    temp=replace(temp,"c3*c4"=>"cc[3,4]")
    for i=1:4
        temp=replace(temp,"c$i"=>"c[$i]")
    end
    for i=10:-1:1
        for j=10:-1:1
                temp=replace(temp,"p$i*q$j"=>"PQ[$i,$j]")
        end
    end
    return temp
end
##
txt=print_term(res[2])
## shape sensetivity?
res=tri_integrate(q1*f1*p1*transpose(gradient(gradient(f1))))*0
for (j,(qq,fq)) in enumerate(zip(q,f))
    for (i,(pp,fp)) in enumerate(zip(p,f))
        println("$j:$i")
        global res
        tmp=fq.subs(z,-z)*transpose(gradient(gradient(fp.subs(z,-z))))*a
        for idx=1:3
            for jdx=1:3
                tmp[idx,jdx]=qq*tmp[idx,jdx].subs(z,0)*pp
            end
        end
        res+=tri_integrate(tmp)
    end
end
##
txt=""
for i =1:3
    for j=1:3
        global txt
txt*="     A2[$i,$j]="
txt*=print_term(res[i,j])
txt*="\n\n"
    end
end
##
tmp=0
test=cc^2*q1*f1*p1*transpose(gradient(f1))
res=tri_integrate(q1*f1*p1*transpose(gradient(f1)))*0
for (j,(qq,fq)) in enumerate(zip(q,f))
    for (i,(pp,fp)) in enumerate(zip(p,f))
        println("$j:$i")
        global res
        tmp=fq*(gradient(fp))*a
        for idx=1:3
                tmp[idx]=qq*tri_integrate(tmp[idx].subs(z,0))*pp
        end
        global tmp
        res+=transpose(tmp)
    end
end
##
txt=""
for i =1:3
        global txt
txt*="     A3[$i]="
txt*=print_term(res[i])
txt*="\n"
end

##
cc^2*q1*f1*p1*gradient(f1)

##Derichlet_Term
res=(cc^2*gradient(f1)*transpose(gradient(f2)))*0
for (j,(qq,fq)) in enumerate(zip(q,f))
    for (i,(pp,fp)) in enumerate(zip(p,f))
                println("$j:$i")
        global res
        tmp=gradient(fq)*transpose(gradient(fp))*a
        for idx=1:3
            for jdx=1:3
                tmp[idx,jdx]=qq*tri_integrate(tmp[idx,jdx].subs(z,0))*pp
            end
        end
        res+=tmp
    end
end
##
txt=""
for i=1:3
    for j=1:3
        global txt
        txt*="M[$i,$j]="*print_term(res[i,j])*"\n\n"
    end
end
