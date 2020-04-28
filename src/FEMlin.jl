import LinearAlgebra
"""mass operator for first order elements"""

function assemble_mass_operator(points,tets)
  #coordinates and number of Points

  #preallocate sparse matrix description
  ii=Array{UInt32}(undef,4^2,length(tets))
  jj=Array{UInt32}(undef,4^2,length(tets))
  mm=Array{ComplexF64}(undef,4^2,length(tets))



  m_unit=[1/60 1/120 1/120 1/120; 1/120 1/60 1/120 1/120; 1/120 1/120 1/60 1/120; 1/120 1/120 1/120 1/60]
  for (idx,smplx) in enumerate(tets)
    tet=points[:,smplx]

    #local coordinate transformation
    J=tet[:,1:3]-tet[:,[4,4,4]]
    #mass_matrix in global coordinates
    m=m_unit*abs(LinearAlgebra.det(J))
    j=repeat(smplx,1,4)
    i=j'
    ii[:,idx]=i[:]
    jj[:,idx]=j[:]
    mm[:,idx]=m[:]
  end

  return ii[:],jj[:],mm[:]

end

"""stiffness operator for first order elements"""
function assemble_stiffness_operator(points, tets, C)
  #preallocate sparse matrix description
  ii=Array{UInt32}(undef,4^2,length(tets))
  jj=Array{UInt32}(undef,4^2,length(tets))
  mm=Array{ComplexF64}(undef,4^2,length(tets))
  m=Array{ComplexF64}(undef,4,4)

  for (idx,(smplx,c)) in enumerate(zip(tets,C))
    tet=points[:,smplx]
    #c=C[idx]
    #local coordinate transformation
    J=tet[:,1:3]-tet[:,[4,4,4]]
    #println(idx,": ",tet)
    J_inv=inv(J)
    A=J_inv*J_inv'

    #stiffness_matrix in global coordinates
    m[1,1]=-A[1,1]*c^2/6
    m[1,2]=-A[1,2]*c^2/6
    m[1,3]=-A[1,3]*c^2/6
    m[1,4]=A[1,1]*c^2/6 + A[1,2]*c^2/6 + A[1,3]*c^2/6
    m[2,1]=m[1,2]
    m[2,2]=-A[2,2]*c^2/6
    m[2,3]=-A[2,3]*c^2/6
    m[2,4]=A[1,2]*c^2/6 + A[2,2]*c^2/6 + A[2,3]*c^2/6
    m[3,1]=m[1,3]
    m[3,2]=m[2,3]
    m[3,3]=-A[3,3]*c^2/6
    m[3,4]=A[1,3]*c^2/6 + A[2,3]*c^2/6 + A[3,3]*c^2/6
    m[4,1]=m[1,4]
    m[4,2]=m[2,4]
    m[4,3]=m[3,4]
    m[4,4]=-A[1,1]*c^2/6 - A[1,2]*c^2/3 - A[1,3]*c^2/3 - A[2,2]*c^2/6 - A[2,3]*c^2/3 - A[3,3]*c^2/6

    m=m*abs(LinearAlgebra.det(J))
    j=repeat(smplx,1,4)
    i=j'
    ii[:,idx]=i[:]
    jj[:,idx]=j[:]
    mm[:,idx]=m[:]
  end

  return ii[:], jj[:], mm[:]
end

"""boundary mass operator for first order elements"""
function assemble_boundary_mass_operator(points, tris, C)
  #preallocate sparse matrix description
  ii=Array{UInt32}(undef,3^2,length(tris))
  jj=Array{UInt32}(undef,3^2,length(tris))
  mm=Array{ComplexF64}(undef,3^2,length(tris))
  m_unit=[1/12 1/24 1/24;
          1/24 1/12 1/12;
          1/24 1/24 1/12]

  for (idx,(smplx,c)) in enumerate(zip(tris,C))
    tri=points[:,smplx]
    #compute triangle surface i.e. determinant of the local coordinate transform
    A=LinearAlgebra.norm(LinearAlgebra.cross(tri[:,1]-tri[:,3],tri[:,2]-tri[:,3]))

    m=m_unit*c*A
    j=repeat(smplx,1,3)#
    i=j'
    ii[:,idx]=i[:]
    jj[:,idx]=j[:]
    mm[:,idx]=m[:]
  end
  return ii[:], jj[:], mm[:]
end

"""volume source vector for first order elements"""
function assemble_volume_source(points,tets)
  #preallocate sparse matrix description
  ii=Array{UInt32}(undef,4,length(tets))
  mm=Array{ComplexF64}(undef,4,length(tets))
  m_unit=[1/24 1/24 1/24 1/24]

  for (idx,smplx) in enumerate(tets)
    tet=points[:,smplx]

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
  ii=smplx
  tet=points[:,smplx]
  J=tet[:,1:3]-tet[:,[4,4,4]]
  J_inv=inv(J)
  grad=[1.0  0.0  0.0;
        0.0  1.0  0.0;
        0.0  0.0  1.0;
       -1.0 -1.0 -1.0]

  mm=grad*J_inv*n_ref
  return ii[:], mm[:]
end

function assemble_boundary_source(points, tris, C,tetmap=nothing)
  #preallocate sparse matrix description
  ii=Array{UInt32}(undef,3,length(tris))
  mm=Array{ComplexF64}(undef,3,length(tris))
  m_unit=[1/6 1/6 1/6]
  for (idx,(smplx,c)) in enumerate(zip(tris,C))
    tri=points[:,smplx]
    #c=C[tetidx] #Todo generate C

    #compute triangle surface i.e. determinant of the local coordinate transform
    A=LinearAlgebra.norm(LinearAlgebra.cross(tri[:,1]-tri[:,3],tri[:,2]-tri[:,3]))

    m=m_unit*c*A #TODO: durch Z
    i=smplx
    ii[:,idx]=i[:]
    mm[:,idx]=m[:]
  end
  return ii[:], mm[:]
end
