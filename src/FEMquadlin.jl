import LinearAlgebra

function assemble_mass_operator(points,tets)
    #preallocate sparse matrix description
    ii=Array{UInt32}(undef,10^2,length(tets))
    jj=Array{UInt32}(undef,10^2,length(tets))
    mm=Array{ComplexF64}(undef,10^2,length(tets))

    m_unit=[1/420 1/2520 1/2520 1/2520 -1/630 -1/630 -1/630 -1/420 -1/420 -1/420 ;
    1/2520 1/420 1/2520 1/2520 -1/630 -1/420 -1/420 -1/630 -1/630 -1/420 ;
    1/2520 1/2520 1/420 1/2520 -1/420 -1/630 -1/420 -1/630 -1/420 -1/630 ;
    1/2520 1/2520 1/2520 1/420 -1/420 -1/420 -1/630 -1/420 -1/630 -1/630 ;
    -1/630 -1/630 -1/420 -1/420 4/315 2/315 2/315 2/315 2/315 1/315 ;
    -1/630 -1/420 -1/630 -1/420 2/315 4/315 2/315 2/315 1/315 2/315 ;
    -1/630 -1/420 -1/420 -1/630 2/315 2/315 4/315 1/315 2/315 2/315 ;
    -1/420 -1/630 -1/630 -1/420 2/315 2/315 1/315 4/315 2/315 2/315 ;
    -1/420 -1/630 -1/420 -1/630 2/315 1/315 2/315 2/315 4/315 2/315 ;
    -1/420 -1/420 -1/630 -1/630 1/315 2/315 2/315 2/315 2/315 4/315]

    for (idx,smplx) in enumerate(tets)
        tet=points[:,smplx[1:4]]
        #local coordinate transformation
        J=tet[:,1:3]-tet[:,[4,4,4]]
        #mass_matrix in global coordinates
        m=m_unit*abs(LinearAlgebra.det(J))
        j=repeat(smplx,1,10)
        i=j'
        ii[:,idx]=i[:]
        jj[:,idx]=j[:]
        mm[:,idx]=m[:]
    end

    return ii[:],jj[:],mm[:]
end

function assemble_stiffness_operator(points,tets,C)
    ii=Array{UInt32}(undef,10^2,length(tets))
    jj=Array{UInt32}(undef,10^2,length(tets))
    mm=Array{ComplexF64}(undef,10^2,length(tets))
    m=Array{ComplexF64}(undef,10,10)

    for (idx,(smplx,c)) in enumerate(zip(tets,C))
      tet=points[:,smplx]
      #local coordinate transformation
      J=tet[:,1:3]-tet[:,[4,4,4]]
      J_inv=inv(J)
      A=J_inv*J_inv'
      #stiffness_matrix in global coordinates
      m[1,1]=-A[1,1]*c^2/10
      m[1,2]=A[1,2]*c^2/30
      m[1,3]=A[1,3]*c^2/30
      m[1,4]=-c^2*(A[1,1] + A[1,2] + A[1,3])/30
      m[1,5]=-4*c^2*(-A[1,1]/120 + A[1,2]/40)
      m[1,6]=-4*c^2*(-A[1,1]/120 + A[1,3]/40)
      m[1,7]=4*c^2*(A[1,1]/30 + A[1,2]/40 + A[1,3]/40)
      m[1,8]=-4*c^2*(-A[1,2]/120 - A[1,3]/120)
      m[1,9]=4*c^2*(-A[1,1]/120 - A[1,3]/120)
      m[1,10]=4*c^2*(-A[1,1]/120 - A[1,2]/120)
      m[2,2]=-A[2,2]*c^2/10
      m[2,3]=A[2,3]*c^2/30
      m[2,4]=-c^2*(A[1,2] + A[2,2] + A[2,3])/30
      m[2,5]=-4*c^2*(A[1,2]/40 - A[2,2]/120)
      m[2,6]=-4*c^2*(-A[1,2]/120 - A[2,3]/120)
      m[2,7]=4*c^2*(-A[2,2]/120 - A[2,3]/120)
      m[2,8]=4*c^2*(A[2,2]/120 - A[2,3]/40)
      m[2,9]=4*c^2*(A[1,2]/40 + A[2,2]/30 + A[2,3]/40)
      m[2,10]=4*c^2*(-A[1,2]/120 - A[2,2]/120)
      m[3,3]=-A[3,3]*c^2/10
      m[3,4]=-c^2*(A[1,3] + A[2,3] + A[3,3])/30
      m[3,5]=-4*c^2*(-A[1,3]/120 - A[2,3]/120)
      m[3,6]=-4*c^2*(A[1,3]/40 - A[3,3]/120)
      m[3,7]=4*c^2*(-A[2,3]/120 - A[3,3]/120)
      m[3,8]=4*c^2*(-A[2,3]/40 + A[3,3]/120)
      m[3,9]=4*c^2*(-A[1,3]/120 - A[3,3]/120)
      m[3,10]=4*c^2*(A[1,3]/40 + A[2,3]/40 + A[3,3]/30)
      m[4,4]=-c^2*(A[1,1] + 2*A[1,2] + 2*A[1,3] + A[2,2] + 2*A[2,3] + A[3,3])/10
      m[4,5]=-4*c^2*(A[1,1]/120 + A[1,2]/60 + A[1,3]/120 + A[2,2]/120 + A[2,3]/120)
      m[4,6]=-4*c^2*(A[1,1]/120 + A[1,2]/120 + A[1,3]/60 + A[2,3]/120 + A[3,3]/120)
      m[4,7]=4*c^2*(A[1,1]/30 + A[1,2]/24 + A[1,3]/24 + A[2,2]/120 + A[2,3]/60 + A[3,3]/120)
      m[4,8]=4*c^2*(-A[1,2]/120 - A[1,3]/120 - A[2,2]/120 - A[2,3]/60 - A[3,3]/120)
      m[4,9]=4*c^2*(A[1,1]/120 + A[1,2]/24 + A[1,3]/60 + A[2,2]/30 + A[2,3]/24 + A[3,3]/120)
      m[4,10]=4*c^2*(A[1,1]/120 + A[1,2]/60 + A[1,3]/24 + A[2,2]/120 + A[2,3]/24 + A[3,3]/30)
      m[5,5]=-16*c^2*(A[1,1]/60 + A[1,2]/60 + A[2,2]/60)
      m[5,6]=-16*c^2*(A[1,1]/120 + A[1,2]/120 + A[1,3]/120 + A[2,3]/60)
      m[5,7]=16*c^2*(A[1,2]/60 + A[1,3]/120 + A[2,2]/60 + A[2,3]/60)
      m[5,8]=8*c^2*(-A[1,2]/60 - A[1,3]/30 - A[2,2]/60 - A[2,3]/60)
      m[5,9]=16*c^2*(A[1,1]/60 + A[1,2]/60 + A[1,3]/60 + A[2,3]/120)
      m[5,10]=16*c^2*(A[1,1]/120 + A[1,2]/60 + A[2,2]/120)
      m[6,6]=-16*c^2*(A[1,1]/60 + A[1,3]/60 + A[3,3]/60)
      m[6,7]=16*c^2*(A[1,2]/120 + A[1,3]/60 + A[2,3]/60 + A[3,3]/60)
      m[6,8]=-8*c^2*(A[1,2]/30 + A[1,3]/60 + A[2,3]/60 + A[3,3]/60)
      m[6,9]=16*c^2*(A[1,1]/120 + A[1,3]/60 + A[3,3]/120)
      m[6,10]=16*c^2*(A[1,1]/60 + A[1,2]/60 + A[1,3]/60 + A[2,3]/120)
      m[7,7]=-16*c^2*(A[1,1]/60 + A[1,2]/60 + A[1,3]/60 + A[2,2]/60 + A[2,3]/30 + A[3,3]/60)
      m[7,8]=8*c^2*(A[2,2]/60 + A[2,3]/30 + A[3,3]/60)
      m[7,9]=-16*c^2*(A[1,2]/60 + A[1,3]/120 + A[2,3]/120 + A[3,3]/120)
      m[7,10]=-16*c^2*(A[1,2]/120 + A[1,3]/60 + A[2,2]/120 + A[2,3]/120)
      m[8,8]=16*c^2*(-A[2,2]/60 - A[2,3]/60 - A[3,3]/60)
      m[8,9]=16*c^2*(A[1,2]/120 + A[1,3]/60 + A[2,3]/60 + A[3,3]/60)
      m[8,10]=16*c^2*(A[1,2]/60 + A[1,3]/120 + A[2,2]/60 + A[2,3]/60)
      m[9,9]=-16*c^2*(A[1,1]/60 + A[1,2]/60 + A[1,3]/30 + A[2,2]/60 + A[2,3]/60 + A[3,3]/60)
      m[9,10]=-16*c^2*(A[1,1]/120 + A[1,2]/120 + A[1,3]/120 + A[2,3]/60)
      m[10,10]=-16*c^2*(A[1,1]/60 + A[1,2]/30 + A[1,3]/60 + A[2,2]/60 + A[2,3]/60 + A[3,3]/60)
      m[2,1]=m[1,2]
      m[3,1]=m[1,3]
      m[3,2]=m[2,3]
      m[4,1]=m[1,4]
      m[4,2]=m[2,4]
      m[4,3]=m[3,4]
      m[5,1]=m[1,5]
      m[5,2]=m[2,5]
      m[5,3]=m[3,5]
      m[5,4]=m[4,5]
      m[6,1]=m[1,6]
      m[6,2]=m[2,6]
      m[6,3]=m[3,6]
      m[6,4]=m[4,6]
      m[6,5]=m[5,6]
      m[7,1]=m[1,7]
      m[7,2]=m[2,7]
      m[7,3]=m[3,7]
      m[7,4]=m[4,7]
      m[7,5]=m[5,7]
      m[7,6]=m[6,7]
      m[8,1]=m[1,8]
      m[8,2]=m[2,8]
      m[8,3]=m[3,8]
      m[8,4]=m[4,8]
      m[8,5]=m[5,8]
      m[8,6]=m[6,8]
      m[8,7]=m[7,8]
      m[9,1]=m[1,9]
      m[9,2]=m[2,9]
      m[9,3]=m[3,9]
      m[9,4]=m[4,9]
      m[9,5]=m[5,9]
      m[9,6]=m[6,9]
      m[9,7]=m[7,9]
      m[9,8]=m[8,9]
      m[10,1]=m[1,10]
      m[10,2]=m[2,10]
      m[10,3]=m[3,10]
      m[10,4]=m[4,10]
      m[10,5]=m[5,10]
      m[10,6]=m[6,10]
      m[10,7]=m[7,10]
      m[10,8]=m[8,10]
      m[10,9]=m[9,10]
      m=m*abs(LinearAlgebra.det(J))
      j=repeat(smplx,1,10)
      i=j'
      ii[:,idx]=i[:]
      jj[:,idx]=j[:]
      mm[:,idx]=m[:]
  end
  return ii[:], jj[:], mm[:]
end

function assemble_boundary_mass_operator(points, tris, C)
  #preallocate sparse matrix description
  ii=Array{UInt32}(undef,6^2,length(tris))
  jj=Array{UInt32}(undef,6^2,length(tris))
  mm=Array{ComplexF64}(undef,6^2,length(tris))
  for (idx,(smplx,c)) in enumerate(zip(tris,C))
    tri=points[:,smplx]
    #compute triangle surface i.e. determinant of the local coordinate transform
    A=LinearAlgebra.norm(LinearAlgebra.cross(tri[:,1]-tri[:,3],tri[:,2]-tri[:,3]))
    m[1,1]=c/60
    m[1,2]=-c/360
    m[1,3]=-c/360
    m[1,4]=0
    m[1,5]=0
    m[1,6]=-c/90
    m[2,2]=c/60
    m[2,3]=-c/360
    m[2,4]=0
    m[2,5]=-c/90
    m[2,6]=0
    m[3,3]=c/60
    m[3,4]=-c/90
    m[3,5]=0
    m[3,6]=0
    m[4,4]=4*c/45
    m[4,5]=2*c/45
    m[4,6]=2*c/45
    m[5,5]=4*c/45
    m[5,6]=2*c/45
    m[6,6]=4*c/45
    m[2,1]=m[1,2]
    m[3,1]=m[1,3]
    m[3,2]=m[2,3]
    m[4,1]=m[1,4]
    m[4,2]=m[2,4]
    m[4,3]=m[3,4]
    m[5,1]=m[1,5]
    m[5,2]=m[2,5]
    m[5,3]=m[3,5]
    m[5,4]=m[4,5]
    m[6,1]=m[1,6]
    m[6,2]=m[2,6]
    m[6,3]=m[3,6]
    m[6,4]=m[4,6]
    m[6,5]=m[5,6]
    m=m*A
    j=repeat(smplx,1,6)#
    i=j'
    ii[:,idx]=i[:]
    jj[:,idx]=j[:]
    mm[:,idx]=m[:]
  end
  return ii[:], jj[:], mm[:]
end

function assemble_volume_source(points,tets)
  #preallocate sparse matrix description
  ii=Array{UInt32}(undef,10,length(tets))
  mm=Array{ComplexF64}(undef,10,length(tets))
  m_unit=[-1/120  -1/120  -1/120  -1/120  1/30  1/30  1/30  1/30  1/30  1/30]

  for (idx,smplx) in enumerate(tets)
    tet=points[:,smplx[1:4]]

    #local coordinate transformation
    J=tet[:,1:3]-tet[:,[4,4,4]]
    #mass_matrix in global coordinates
    m=m_unit*abs(LinearAlgebra.det(J))
    ii[:,idx]=smplx[:]
    mm[:,idx]=m[:]
  end

  return ii[:], mm[:]
end

function assemble_gradient_source(points,smplx, x_ref, n_ref)
    tet=points[:,smplx[1:4]]
    ii=smplx
    #local coordinate transformation
    J=tet[:,1:3]-tet[:,[4,4,4]]
    J_inv=inv(J)
    x,y,z=J_inv*(x_ref-points[:,smplx[4]])
    grad=[4*x-1 0 0 ;
          0 4*y-1 0 ;
          0 0 4*z-1 ;
          4*x+4*y+4*z-3 4*x+4*y+4*z-3 4*x+4*y+4*z-3 ;
          4*y 4*x 0 ;
          4*z 0 4*x ;
          -8*x-4*y-4*z+4 -4*x -4*x ;
          0 4*z 4*y ;
          -4*y -4*x-8*y-4*z+4 -4*y ;
          -4*z -4*z -4*x-4*y-8*z+4 ]
    mm=grad*J_inv*n_ref
    return ii[:], mm[:]
end

function assemble_boundary_source(points, tris, C,tetmap=nothing)
  #preallocate sparse matrix description
  ii=Array{UInt32}(undef,6,length(tris))
  mm=Array{ComplexF64}(undef,6,length(tris))
  m_unit=[ 0 0 0 1/6 1/6 1/6]
  for (idx,(smplx,c)) in enumerate(zip(tris,C))
    tri=points[:,smplx]
    #compute triangle surface i.e. determinant of the local coordinate transform
    A=LinearAlgebra.norm(LinearAlgebra.cross(tri[:,1]-tri[:,3],tri[:,2]-tri[:,3]))

    m=m_unit*c*A #TODO: durch Z
    i=smplx
    ii[:,idx]=i[:]
    mm[:,idx]=m[:]
  end
  return ii[:], mm[:]
end
