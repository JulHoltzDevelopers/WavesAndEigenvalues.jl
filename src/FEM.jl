#local coordinate transformation type
struct CooTrafo
        trafo
        inv
        det
end
#constructor
function CooTrafo(X)
        dim=size(X)
        J=Array{Float64}(undef,dim[1],dim[1])
        J[:,1:dim[2]-1]=X[:,1:end-1].-X[:,end]
        if dim[2]==3 #surface triangle #TODO:implement lines and all dimensions
                n=LinearAlgebra.cross(J[:,1],J[:,2])
                J[:,end]=n./LinearAlgebra.norm(n)
        end
        #println("###HERE###")
        #println(X)
        #println(J)
        Jinv=LinearAlgebra.inv(J)
        det=LinearAlgebra.det(J)
        CooTrafo(J,Jinv,det)
end
function create_indices(smplx)
    ii=Array{Int64}(undef,length(smplx),length(smplx))
    for i=1:length(smplx)
        for j=1:length(smplx)
            ii[i,j]=smplx[i]
        end
    end

    jj=ii'
    return ii, jj
end

# Functions to compute local finite element matrices on simplices. The function
# names follow the pattern `sAB[d]v[d]uC[[d]cD]`. Where
# - `A`: is the number of vertices of the simplex, e.g. 4 for a tetrahedron
# - `B`: is the number of space dimension
# - `C`: the order of the ansatz functions
# - `D`: the interpolation order of the (optional) coefficient function
#
# optional `d` indicate a partial derivative.
#
# Example: Lets assume you want to discretize the term ∂u(x,y,z)/∂x*c(x,y,z)
# Then the local weak form gives rise to the integral
# ∫∫∫_(Tet) v*∂u(x,y,z)/∂x*c(x,y,z) dxdydz  where `Tet` is the tetrahedron.
# Furthermore, trial and test functions u and v should be of second order
# while the coefficient function should be interpolated linearly.
# Then, the function that returns the local discretization matrix of this weak
# form is s43vdu2c1. (The direction of the partial derivative is specified in
# the function call.)

#TODO: multiple dispatch instead of/additionally to function names

## mass matrices
function s33v1u1(J::CooTrafo)
     M=[1/12 1/24 1/24;
        1/24 1/12 1/24;
        1/24 1/24 1/12]

    return M*abs(J.det)
end


function s43v1u1(J::CooTrafo)
    M = [1/60 1/120 1/120 1/120;
          1/120 1/60 1/120 1/120;
          1/120 1/120 1/60 1/120;
          1/120 1/120 1/120 1/60]
    return M*abs(J.det)
end

function s43v1u1c1(J::CooTrafo,c)
        c1,c2,c3,c4=c
        M=Array{ComplexF64}(undef,4,4)
        M[1,1]=c1/120 + c2/360 + c3/360 + c4/360
        M[1,2]=c1/360 + c2/360 + c3/720 + c4/720
        M[1,3]=c1/360 + c2/720 + c3/360 + c4/720
        M[1,4]=c1/360 + c2/720 + c3/720 + c4/360
        M[2,2]=c1/360 + c2/120 + c3/360 + c4/360
        M[2,3]=c1/720 + c2/360 + c3/360 + c4/720
        M[2,4]=c1/720 + c2/360 + c3/720 + c4/360
        M[3,3]=c1/360 + c2/360 + c3/120 + c4/360
        M[3,4]=c1/720 + c2/720 + c3/360 + c4/360
        M[4,4]=c1/360 + c2/360 + c3/360 + c4/120
        M[2,1]=M[1,2]
        M[3,1]=M[1,3]
        M[4,1]=M[1,4]
        M[3,2]=M[2,3]
        M[4,2]=M[2,4]
        M[4,3]=M[3,4]
        return M*abs(J.det)
end

function s43v2u2(J::CooTrafo)
    M = [1/420 1/2520 1/2520 1/2520 -1/630 -1/630 -1/630 -1/420 -1/420 -1/420;
            1/2520 1/420 1/2520 1/2520 -1/630 -1/420 -1/420 -1/630 -1/630 -1/420;
            1/2520 1/2520 1/420 1/2520 -1/420 -1/630 -1/420 -1/630 -1/420 -1/630;
            1/2520 1/2520 1/2520 1/420 -1/420 -1/420 -1/630 -1/420 -1/630 -1/630;
            -1/630 -1/630 -1/420 -1/420 4/315 2/315 2/315 2/315 2/315 1/315;
            -1/630 -1/420 -1/630 -1/420 2/315 4/315 2/315 2/315 1/315 2/315;
            -1/630 -1/420 -1/420 -1/630 2/315 2/315 4/315 1/315 2/315 2/315;
            -1/420 -1/630 -1/630 -1/420 2/315 2/315 1/315 4/315 2/315 2/315;
            -1/420 -1/630 -1/420 -1/630 2/315 1/315 2/315 2/315 4/315 2/315;
            -1/420 -1/420 -1/630 -1/630 1/315 2/315 2/315 2/315 2/315 4/315]
    return M*abs(J.det)
end

## partial derivatives
function s43v1du1(J,d)
        M1=     [1/24 0 0 -1/24;
                1/24 0 0 -1/24;
                1/24 0 0 -1/24;
                1/24 0 0 -1/24]

        M2=     [0 1/24 0 -1/24;
                0 1/24 0 -1/24;
                0 1/24 0 -1/24;
                0 1/24 0 -1/24]

        M3=     [0 0 1/24 -1/24;
                0 0 1/24 -1/24;
                0 0 1/24 -1/24;
                0 0 1/24 -1/24]
        return (M1.*J.inv[1,d].+M2.*J.inv[2,d].+M3.*J.inv[3,d])*abs(J.det)
end

s43dv1u1(J,d)=s43v1du1(J,d)'

function s43v1du1c1(J::CooTrafo,c,d)
        c1,c2,c3,c4=c
        M1=Array{ComplexF64}(undef,4,4)
        M2=Array{ComplexF64}(undef,4,4)
        M3=Array{ComplexF64}(undef,4,4)

        M1[1,1]=c1/60 + c2/120 + c3/120 + c4/120
        M1[1,2]=0
        M1[1,3]=0
        M1[1,4]=-c1/60 - c2/120 - c3/120 - c4/120
        M1[2,1]=c1/120 + c2/60 + c3/120 + c4/120
        M1[2,2]=0
        M1[2,3]=0
        M1[2,4]=-c1/120 - c2/60 - c3/120 - c4/120
        M1[3,1]=c1/120 + c2/120 + c3/60 + c4/120
        M1[3,2]=0
        M1[3,3]=0
        M1[3,4]=-c1/120 - c2/120 - c3/60 - c4/120
        M1[4,1]=c1/120 + c2/120 + c3/120 + c4/60
        M1[4,2]=0
        M1[4,3]=0
        M1[4,4]=-c1/120 - c2/120 - c3/120 - c4/60

        M2[1,1]=0
        M2[1,2]=c1/60 + c2/120 + c3/120 + c4/120
        M2[1,3]=0
        M2[1,4]=-c1/60 - c2/120 - c3/120 - c4/120
        M2[2,1]=0
        M2[2,2]=c1/120 + c2/60 + c3/120 + c4/120
        M2[2,3]=0
        M2[2,4]=-c1/120 - c2/60 - c3/120 - c4/120
        M2[3,1]=0
        M2[3,2]=c1/120 + c2/120 + c3/60 + c4/120
        M2[3,3]=0
        M2[3,4]=-c1/120 - c2/120 - c3/60 - c4/120
        M2[4,1]=0
        M2[4,2]=c1/120 + c2/120 + c3/120 + c4/60
        M2[4,3]=0
        M2[4,4]=-c1/120 - c2/120 - c3/120 - c4/60

        M3[1,1]=0
        M3[1,2]=0
        M3[1,3]=c1/60 + c2/120 + c3/120 + c4/120
        M3[1,4]=-c1/60 - c2/120 - c3/120 - c4/120
        M3[2,1]=0
        M3[2,2]=0
        M3[2,3]=c1/120 + c2/60 + c3/120 + c4/120
        M3[2,4]=-c1/120 - c2/60 - c3/120 - c4/120
        M3[3,1]=0
        M3[3,2]=0
        M3[3,3]=c1/120 + c2/120 + c3/60 + c4/120
        M3[3,4]=-c1/120 - c2/120 - c3/60 - c4/120
        M3[4,1]=0
        M3[4,2]=0
        M3[4,3]=c1/120 + c2/120 + c3/120 + c4/60
        M3[4,4]=-c1/120 - c2/120 - c3/120 - c4/60

        return (M1.*J.inv[1,d].+M2.*J.inv[2,d].+M3.*J.inv[3,d])*abs(J.det)
end

s43dv1u1c1(J::CooTrafo,c,d)=s43v1du1c1(J::CooTrafo,c,d)'


function s43v1u1dc1(J::CooTrafo,c,d)
        c1,c2,c3,c4=c
        M = [1/60 1/120 1/120 1/120;
              1/120 1/60 1/120 1/120;
              1/120 1/120 1/60 1/120;
              1/120 1/120 1/120 1/60]

        return ([c1-c4, c2-c4, c3-c4]*J.inv[:,d]).*M*abs(J.det)
end

# nabla operations
function s43nv1nu1(J::CooTrafo)
        A=J.inv*J.inv'
        M=Array{Float64}(undef,4,4)
        M[1,1]=A[1,1]/6
        M[1,2]=A[1,2]/6
        M[1,3]=A[1,3]/6
        M[1,4]=-A[1,1]/6 - A[1,2]/6 - A[1,3]/6
        M[2,2]=A[2,2]/6
        M[2,3]=A[2,3]/6
        M[2,4]=-A[1,2]/6 - A[2,2]/6 - A[2,3]/6
        M[3,3]=A[3,3]/6
        M[3,4]=-A[1,3]/6 - A[2,3]/6 - A[3,3]/6
        M[4,4]=A[1,1]/6 + A[1,2]/3 + A[1,3]/3 + A[2,2]/6 + A[2,3]/3 + A[3,3]/6
        M[2,1]=M[1,2]
        M[3,1]=M[1,3]
        M[4,1]=M[1,4]
        M[3,2]=M[2,3]
        M[4,2]=M[2,4]
        M[4,3]=M[3,4]
        return M*abs(J.det)
end

## source vectors

function s33v1(J)
        return [1/6 1/6 1/6]*abs(J.det)
end
