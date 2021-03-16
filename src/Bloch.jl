export bloch_expand


function blochify(ii,jj,mm, naxis,nxbloch,nsector,naxis_ln,nsector_ln,N_points; axis = true)
    #nsector=naxis+nxbloch+nbody+nxsymmetry+nbody
    blochshift=nsector-naxis
    blochshift_ln=nsector_ln-naxis_ln# blochshift is necessary because point dofs are reduced by blochshift
    MM=ComplexF64[]
    MM_plus=ComplexF64[]
    MM_minus=ComplexF64[]
    MM_plus_axis=ComplexF64[]
    MM_minus_axis=ComplexF64[]
    MM_axis=ComplexF64[]

    II=UInt32[]
    JJ=UInt32[]
    II_plus=UInt32[]#TODO: consider preallocation
    JJ_plus=UInt32[]
    II_minus=UInt32[]
    JJ_minus=UInt32[]
    II_axis=UInt32[]
    JJ_axis=UInt32[]
    II_plus_axis=UInt32[]
    JJ_plus_axis=UInt32[]
    II_minus_axis=UInt32[]
    JJ_minus_axis=UInt32[]
    #println("$(length(ii)) $(length(jj)) $(length(mm))")
    for (i,j,m) in zip(ii,jj,mm)
        if i<=N_points #deal with point dof
            i_check = i>nsector
            # map Bloch_ref to Bloch_img
            if i_check
                i-=blochshift
            end
        else #deal with line dof
            i_check = i> nsector_ln
            if i_check
                i-=blochshift_ln
            end
        end
        if j<=N_points #deal with point dof
            j_check = j>nsector
            # map Bloch_ref to Bloch_img
            if j_check
                j-=blochshift
            end
        else #deal with line dof
            j_check = j> nsector_ln
            if j_check
                j-=blochshift_ln
            end
        end

        if axis && ( i<=naxis || j<=naxis || N_points<i<=naxis_ln || N_points<j<=naxis_ln) #&& !(i<=naxis && j<=naxis) && !(N_points>i<=naxis_ln && N_points>j<=naxis_ln) &&!(i<=naxis && N_points>j<=naxis_ln) && !(N_points>i<=naxis_ln && j<=naxis)
            axis_check=true
        else
            axis_check=false
        end

        #account for reduced number of point dofs
        if i>N_points
            i-=nxbloch
        end
        if j>N_points
            j-=nxbloch
        end

        #sort into operators matrices
        if (!i_check && !j_check) || (i_check && j_check)  #no manipulation of matrix entries
            if axis_check
                append!(II_axis,i)
                append!(JJ_axis,j)
                append!(MM_axis,m)
            else
                append!(II,i)
                append!(JJ,j)
                append!(MM,m)
            end

        elseif !i_check && j_check
            if axis_check
                append!(II_plus_axis,i)
                append!(JJ_plus_axis,j)
                append!(MM_plus_axis,m)
            else
                append!(II_plus,i)
                append!(JJ_plus,j)
                append!(MM_plus,m)
            end

        elseif i_check && !j_check
            if axis_check
                append!(II_minus_axis,i)
                append!(JJ_minus_axis,j)
                append!(MM_minus_axis,m)
            else
                append!(II_minus,i)
                append!(JJ_minus,j)
                append!(MM_minus,m)
            end
        else
            println("ERROR: in blochification")
            return nothing
        end
    end

    if naxis==0
        return (II,II_plus,II_minus), (JJ,JJ_plus,JJ_minus), (MM, MM_plus, MM_minus)
    else
        return (II, II_plus, II_minus, II_axis, II_plus_axis, II_minus_axis), (JJ, JJ_plus, JJ_minus, JJ_axis, JJ_plus_axis, JJ_minus_axis), (MM, MM_plus, MM_minus, MM_axis, MM_plus_axis, MM_minus_axis)
    end
end

"""
    v=bloch_expand(mesh::Mesh,sol::Solution,b=:b)

Expand solution vector `sol.v` on mesh `mesh`with Bloch wave number
`sol.params[b]` and return it as `v`.
"""
function bloch_expand(mesh::Mesh,sol::Solution,b=:b)
    naxis=mesh.dos.naxis
    nxsector=mesh.dos.nxsector
    DOS=mesh.dos.DOS
    v=zeros(ComplexF64,naxis+nxsector*DOS)
    v[1:naxis]=sol.v[1:naxis]
    B=sol.params[b]
    for s=0:DOS-1
        v[naxis+1+s*nxsector:naxis+(s+1)*nxsector]=sol.v[naxis+1:naxis+nxsector].*exp(+2.0im*pi/DOS*B*s)
    end
    return v
end

function bloch_expand(mesh::Mesh,vec::Array,b::Real=0)
    naxis=mesh.dos.naxis
    nxsector=mesh.dos.nxsector
    DOS=mesh.dos.DOS
    v=zeros(ComplexF64,naxis+nxsector*DOS)
    v[1:naxis]=vec[1:naxis]
    for s=0:DOS-1
        v[naxis+1+s*nxsector:naxis+(s+1)*nxsector]=vec[naxis+1:naxis+nxsector].*exp(+2.0im*pi/DOS*b*s)
    end
    return v
end
