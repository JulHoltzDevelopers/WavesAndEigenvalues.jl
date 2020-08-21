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
#barycentric coordinate vector
X=[x,y,z,a]
# shape function list
F= [[x,y,z,a], [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10 ],
    [fh01, fh02,fh03,fh04,fh05,fh06,fh07,fh08,fh09,fh10,fh11,fh12,fh13,fh14,fh15,fh16,fh17,fh18,fh19,fh20],
    [fh01, fh02,fh03,fh04,fr05,fr06,fr07,fr08,fr09,fr10,fr11,fr12,fr13,fr14,fr15,fr16,fh17,fh18,fh19,fh20]]

#coefficient lists
C= [[c1*x, c2*y, c3*z, c4*a,],[c1*f1, c2*f2, c3*f3, c4*f4, c5*f5, c6*f6, c7*f7, c8*f8, c9*f9, c10*f10 ],]
b= [[b1*x, b2*y, b3*z, b4*a,],[b1*f1, b2*f2, b3*f3, b4*f4, b5*f5, b6*f6, b7*f7, b8*f8, b9*f9, b10*f10 ],]

#f = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10 ]
#ff=[f1,f4,f5,f7]

#surface point list
surf=[[1,2,4],[1,2,4,5,7,9],[1 2 4 5 6 8 9 10 12 13 14 16 19],[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]]


##
# function check_hermitian_mat(F)
#     M=Array{Any}(undef,20,20)
#     for (i,f) in enumerate(F)
#         for (j,pnt) in enumerate([(1,0,0), (0,1,0), (0,0,1), (0,0,0), (0,1/3,1/3), (1/3,0,1/3), (1/3,1/3,0), (1/3,1/3,1/3)])
#             xx,yy,zz=pnt
#             res=f.subs(x,xx).subs(y,yy).subs(z,zz).simplify()
#             #println("$pnt ---> $res")
#
#             res=diff(f,x).subs(x,xx).subs(y,yy).subs(z,zz).simplify()
#             if j >4
#                 idx=j+12
#             else
#                 idx=j
#             end
#             M[idx,i]=res
#
#             if j<=4
#                 res=diff(f,x).subs(x,xx).subs(y,yy).subs(z,zz).simplify()
#                 #println("dx $pnt ---> $res")
#                 M[j+4,i]=res
#                 res=diff(f,y).subs(x,xx).subs(y,yy).subs(z,zz).simplify()
#                 #println("dy $pnt ---> $res")
#                 M[j+8,i]=res
#                 res=diff(f,z).subs(x,xx).subs(y,yy).subs(z,zz).simplify()
#                 #println("dz $pnt ---> $res")
#                 M[j+12,i]=res
#             end
#         end
#     end
#     return M
# end
# TT=check_hermitian_mat(F[3])
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
# function gradient(f)
#     N=length(f)
#     return transpose(SymPy.Matrix([SymPy.diff(f,x) SymPy.diff(f,y) SymPy.diff(f,z)]))
# end

# function build_gradient_matrix(f)
#     M=gradient(f[0])
#     for i = 1:N-1
#         M=M.col_insert(i,gradient(f[i]).T)
#     end
#     return M
# end

function tet_integrate(f)
    if size(f)!=()
        f=f[1]
    end
    return SymPy.integrate(SymPy.integrate(SymPy.integrate(f, (x, 0, 1 - z -y)), (y, 0, 1 - z)), (z, 0, 1))
end

function tri_integrate(f)
    return SymPy.integrate(SymPy.integrate(f, (x, 0, 1 - y)), (y, 0, 1))
end
