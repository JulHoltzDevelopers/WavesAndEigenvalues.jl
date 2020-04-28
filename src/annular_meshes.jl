## geometry functions
"""
    pln=three_points_to_plane(A)

Compute equation for a plane from the three points defined in A. The plane is parameterized
as `a*x+b*y+c*z+d==0` with the coefficients returned as `pln=[a,b,c,d]`.

#Arguments
- A::3×3-Array : Array containing the three points defining the plane as columns.

#Notes
Code adapted from https://www.geeksforgeeks.org/program-to-find-equation-of-a-plane-passing-through-3-points/
"""
function three_points_to_plane(A)
    A,B,C=A[:,1],A[:,2],A[:,3]
    u=A-C
    v=B-C
    a=LinearAlgebra.cross(u,v)
    a/=LinearAlgebra.norm(a) #normalize
    d=-LinearAlgebra.dot(a,C)
    if d<1E-7 #TODO: fidn a better fix for the intersection problem
        d=0.0
    end
    a,b,c=a
    return [a,b,c,d]
end


"""
    check = is_point_in_plane(pnt, pln)

Check whether point "pnt" is in plane . Note that there is currently no handling of round off errors.
"""
function is_point_in_plane(pnt,pln)
    #TODO: implement round off handling
    return pln[1]*pnt[1]+pln[2]*pnt[2]+pln[3]*pnt[3]+pln[4]
end
"""
    p=reflect_point_at_plane(pnt,pln)

Reflect point `pnt` at plane `pln` and return as `p`.

#Notes
Code adapted from https://www.geeksforgeeks.org/mirror-of-a-point-through-a-3-d-plane/
"""
function reflect_point_at_plane(pnt,pln)
    a,d=pln[1:3],pln[4]
    k=(-LinearAlgebra.dot(a,pnt)-d)#/(LinearAlgebra.dot(a,a))
    p=a.*k.+pnt #foot of the perpendicular
    p=2*p-pnt #reflected point
    return p
end
"""
    foot=find_foot_of_perpendicular(pnt,pln)

Compute the position of the foot of the perpendicular from the point `pnt` to the plane `pln`.

#Notes
Code adapted from https://www.geeksforgeeks.org/mirror-of-a-point-through-a-3-d-plane/
"""
function find_foot_of_perpendicular(pnt,pln)
    a,d=pln[1:3],pln[4]
    k=(-LinearAlgebra.dot(a,pnt)-d)#/(LinearAlgebra.dot(a,a))
    foot=a.*k.+pnt #foot of the perpendicular
    return foot
end

"""
 make_normal_outwards!(pln,testpoint)

Ensure that the parameterization of the plane "pln" has a normal that is directed twoards `testpoint`.
If necessary reparametrize `pln`.
"""
function make_normal_outwards!(pln,testpoint)
    foot=find_foot_of_perpendicular(testpoint,pln)
    scalar_product=LinearAlgebra.dot(pln[1:3],testpoint-foot)
    pln*=-sign(scalar_product)
    return
end

"""
    p,n = find_intersection_of_two_planes(pln1,pln2)

Find a parametrization of the axis defined by the intersection of the planes `pln1`and `pln2`.
The axis is parameterized by a point `p` and direction `n`.

#Notes
The algorithm follows John Krumm's solution for finding a common point on two intersecting planes as it is
explained in https://math.stackexchange.com/questions/475953/how-to-calculate-the-intersection-of-two-planes
To simplify the code, here the point p0 is the origin.
"""
function find_intersection_of_two_planes(pln1,pln2)
    n=LinearAlgebra.cross(pln1[1:3],pln2[1:3]) #axis vector
    n/=LinearAlgebra.norm(n)
    #TODO: check whether vector is nonzero and planes truly intersect!
    #generate points on the planes closest to origin.
    p1=find_foot_of_perpendicular([0.,0.,0.],pln1)
    p2=find_foot_of_perpendicular([0.,0.,0.],pln2)
    #Linear system from Lagrange optimization
    M=[2.0 0.0 0.0 pln1[1] pln2[1];
       2.0 0.0 0.0 pln1[2] pln2[2];
       2.0 0.0 0.0 pln1[3] pln2[3];
     pln1[1] pln1[2] pln1[3] 0 0;
     pln2[1] pln2[2] pln2[3] 0 0]

     rhs=[0;0;0; LinearAlgebra.dot(p1,pln1[1:3]); LinearAlgebra.dot(p2,pln2[1:3])]
     p=rhs\M
     return p[1:3],n
end

"""
    R = create_rotation_matrix_around_axis(n,α)

Compute the rotation matrix around the directional vector `n`rotating by an angle `α` (in rad).
"""
function create_rotation_matrix_around_axis(n,α)
    R=[n[1]^2*(1-cos(α))+cos(α) n[1]*n[2]*(1-cos(α))-n[3]*sin(α) n[1]*n[3]*(1-cos(α))+n[2]*sin(α);
       n[2]*n[1]*(1-cos(α))+n[3]*sin(α) n[2]^2*(1-cos(α))+cos(α) n[2]*n[3]*(1-cos(α))-n[1]*sin(α);
       n[3]*n[1]*(1-cos(α))-n[2]*sin(α) n[3]*n[2]*(1-cos(α))+n[1]*sin(α) n[3]^2*(1-cos(α))+cos(α)]
       return R
end

## mesh book keeping
function find_testpoint_idx(tri, tet)
    testpoint_idx=0
    for idx in tet
        if !(idx in tri)
            testpoint_idx=idx
        end
    end
    return testpoint_idx
end

##

"""
    full_mesh=create_full_mesh(mesh::Mesh, doms; sym_name="Symmetry", blch_name="Bloch", unit=false)

Create full mesh from half cell represented in `mesh`.

# Arguments
- `mesh::Mesh`: mesh representing the half cell. The mesh must span a sector of `2π/2N`, where `N` is the (integer) degree of symmetry of the full mesh.
- `doms::List`: list of 2-tuples containing the domain names and their `copy_degree` (see notes below).
- `sym_name::String="Symmetry"`: name of the domain of `mesh` that forms the symmetry plane of the half-mesh.
- `blch_name::String="Bloch"`: name of the domain of `mesh` that forms the remaining azimuthal plane of the half-mesh.
- `unit::Bool=false`: toggle whether extend the mesh to unit cell only.

# Returns
- `full_mesh::Mesh`: representation of the full mesh.

# Notes
The routine copies only domains that are specified  in `doms`. These domains are extended according to the specified `copy_degree`.
The following are available:

- `:full`: extent the domain and save it under the same name in `full_mesh`.
- `:unit`: extent the domain and save the individual unit cells labeled from `0` to `N-1` in `full_mesh`.#
- `:half`: extent the domain and save the individual unit cells as half cells labeled from `0` to `N-1` in full_mesh where one half-cell  contains `_img` in its name.

Multiple styles can be mixed in one domain specification. An example for `doms` would be
    doms=[("Interior", :full), ("Outlet", :unit), ("Flame", :half)]
"""
function create_full_mesh(mesh::Mesh, doms; sym_name="Symmetry", blch_name="Bloch", unit=false)
    npoints=size(mesh.points,2)
    # axis=[]
    # for smplx_idx in mesh.domains["Centeraxis"]["simplices"]
    #     for pnt in mesh.lines[smplx_idx]
    #         insert_smplx!(axis,pnt)
    #     end
    # end

    bloch=[]
    for smplx_idx in mesh.domains[blch_name]["simplices"]
        for pnt in mesh.triangles[smplx_idx]
            insert_smplx!(bloch,pnt)
        end
    end

    symmetry=[]
    for smplx_idx in mesh.domains[sym_name]["simplices"]
        for pnt in mesh.triangles[smplx_idx]
            insert_smplx!(symmetry,pnt)
        end
    end

    axis=intersect(bloch,symmetry)

    #check whether two methods yield same result
    #caxis==axis
    #resort points
    new_order=[]

    for idx in axis
        append!(new_order,idx)
    end

    for idx in bloch
        if !(idx in new_order)
            append!(new_order,idx)
        end
    end

    for idx in 1:npoints
        if !(idx in new_order) && !(idx in symmetry)
            append!(new_order,idx)
        end
    end

    for idx in symmetry
        if !(idx in new_order)
            append!(new_order,idx)
        end
    end

    trace_order=zeros(Int64,npoints) #TODO: loop is not so nice but does the job...
    for idx = 1: length(new_order)
        trace_order[idx]=findfirst(x->x==idx,new_order)
    end
    #next mirror at symmetry

    ##
    # function link_triangles_to_tetrahedra(mesh::Mesh)
    #     tetmap=Array{UInt32,1}(undef,length(mesh.triangles))
    #     for (idx,tri) in enumerate(mesh.triangles)
    #         for (idt::UInt32,tet) in enumerate(mesh.tetrahedra)
    #             if all([i in tet for i=tri])
    #                 tetmap[idx]=idt
    #                 break
    #             end
    #         end
    #     end
    #     return tetmap
    # end




    # compute symmetry plane
    link_triangles_to_tetrahedra!(mesh)
    tetmap=mesh.tri2tet

    smplx_idx=mesh.domains[sym_name]["simplices"][1]
    smplx=mesh.triangles[smplx_idx]
    tri=mesh.points[:,smplx]
    pln=three_points_to_plane(tri)
    pnt=find_foot_of_perpendicular([0.,0.,0.],pln)
    testpoint_idx=find_testpoint_idx(smplx,mesh.tetrahedra[tetmap[smplx_idx]])
    make_normal_outwards!(pln,mesh.points[:,testpoint_idx])


    # compute bloch plane
    smplx_idx=mesh.domains[blch_name]["simplices"][2]
    smplx=mesh.triangles[smplx_idx]
    tri2=mesh.points[:,smplx]
    bpln=three_points_to_plane(tri2)
    testpoint_idx=find_testpoint_idx(smplx,mesh.tetrahedra[tetmap[smplx_idx]])
    make_normal_outwards!(bpln,mesh.points[:,testpoint_idx])
    #little tests
    #check=is_point_in_plane(tri2[:,2],pln)
    #tri_refl=reflect_point_at_plane(tri2[:,2],pln)


    # do the reflection
    nbloch=length(bloch)
    naxis=length(axis)
    nsymmetry=length(symmetry)
    nxsymmetry=nsymmetry-naxis
    nxbloch=nbloch-naxis
    nbody=npoints-nbloch-nxsymmetry
    shiftbody=npoints-nbloch
    shiftbloch=nxbloch+nbody+nxsymmetry+nbody
    #initialize point array
    points=zeros(3,2*npoints-nsymmetry)
    points[:,1:npoints]=mesh.points[:,new_order]
    # reflect body points
    for idx=nbloch+1:npoints-nxsymmetry
        pnt=points[:,idx]
        points[:,idx+shiftbody]=reflect_point_at_plane(pnt,pln)
    end
    #reflect exclusive bloch points
    for idx=naxis+1:nbloch
        pnt=points[:,idx]
        points[:,idx+shiftbloch]=reflect_point_at_plane(pnt,pln)
    end

    #get degree of symmetry
    phi=LinearAlgebra.dot(pln[1:3],-bpln[1:3])#/(LinearAlgebra.norm(pln[1:3])*LinearAlgebra.norm(bpln[1:3]))
    phi=acos(phi)##TODO: get direction
    DOS=round(Int,pi/phi)
    #create/virtualize full mesh

    p,n=find_intersection_of_two_planes(pln,bpln)
    phi=2pi/DOS
    R=create_rotation_matrix_around_axis(n,phi)
    #
    #rot_points=R*(points.-p).+p #Test
    #Full mesh
    nfpoints=naxis+(nxbloch+nbody+nxsymmetry+nbody)*DOS

    if unit
        fpoints=points
        DOS_lim=1
        file_txt="unit from $(mesh.file)"
    else
        DOS_lim=DOS
        fpoints=zeros(3,nfpoints)

        fpoints[:,1:naxis+nxbloch+nbody+nxsymmetry+nbody]=points[:,1:naxis+nxbloch+nbody+nxsymmetry+nbody]
        for idx in 1:DOS-1
            R=create_rotation_matrix_around_axis(n,idx*phi)
            fpoints[:,naxis+(nxbloch+nbody+nxsymmetry+nbody)*idx+1:naxis+(nxbloch+nbody+nxsymmetry+nbody)*(idx+1)]=R*(points[:,naxis+1:naxis+(nxbloch+nbody+nxsymmetry+nbody)].-p).+p
        end
        file_txt="extended from $(mesh.file)"
    end
    domains=Dict()
    lines=Array{UInt32,1}[]
    triangles=Array{UInt32,1}[]


    function get_reflected_index(idx)
        if idx <= naxis
            idx=idx
        elseif idx <= naxis+nxbloch
            idx=idx+shiftbloch
        elseif idx <= naxis+nxbloch+nbody
                idx=idx+shiftbody
        elseif idx <=  naxis+nxbloch+nbody+nxsymmetry
            idx=idx
        else
            idx=0 #error flag
        end
        return idx
    end


    shiftrot=nxbloch+nbody+nxsymmetry+nbody
    function get_rotated_index(idx,sector)
        if idx <= naxis
            idx=idx
        else
            idx=idx+shiftrot*sector #sector counting starts from zero
        end

        #modulo
        if idx>naxis+shiftrot*DOS
            idx=idx-shiftrot*DOS
        end
        return idx
    end

    tetrahedra=Array{UInt32,1}[]


    for (smplx_idx,tet) in enumerate(mesh.tetrahedra)
        #println(idx," ", length(tetrahedra))
        tet=trace_order[tet]
        rtet=get_reflected_index.(tet)

        for sector=0:DOS_lim-1
            insert_smplx!(tetrahedra,get_rotated_index.(tet, sector))
            insert_smplx!(tetrahedra,get_rotated_index.(rtet, sector))
        end
    end


    triangles=Array{UInt32,1}[]
    for (smplx_idx,smplx) in enumerate(mesh.triangles)
        if !(smplx_idx in mesh.domains[sym_name]["simplices"]) && (!(smplx_idx in mesh.domains[blch_name]["simplices"]) || unit)
            smplx=trace_order[smplx]
            rsmplx=get_reflected_index.(smplx)
            for sector=0:DOS_lim-1
                insert_smplx!(triangles,get_rotated_index.(smplx, sector))
                insert_smplx!(triangles,get_rotated_index.(rsmplx, sector))
            end
        end
    end

    #build domains
    domains=Dict()
    #build unit cell

    #build full circumference #TODO: virtualize domains to save memeory and building time
    for (dom, copy_degree) in doms
        dim=mesh.domains[dom]["dimension"]
        if dim == 3
            list_of_simplices=mesh.tetrahedra
            full_list_of_simplices=tetrahedra
        elseif dim == 2
            list_of_simplices=mesh.triangles
            full_list_of_simplices=triangles
        end

        if copy_degree==:full
            domains[dom]=Dict()
            domains[dom]["dimension"]=dim
            domains[dom]["simplices"]=UInt64[]
        elseif copy_degree==:unit
            for sector=0:DOS_lim-1
                domains[dom*"#$sector"]=Dict()
                domains[dom*"#$sector"]["dimension"]=dim
                domains[dom*"#$sector"]["simplices"]=UInt64[]
            end
        elseif copy_degree==:half
            for sector=0:DOS_lim-1
                domains[dom*"#$sector.0"]=Dict()
                domains[dom*"#$sector.0"]["dimension"]=dim
                domains[dom*"#$sector.0"]["simplices"]=UInt64[]
                domains[dom*"#$sector.1"]=Dict()
                domains[dom*"#$sector.1"]["dimension"]=dim
                domains[dom*"#$sector.1"]["simplices"]=UInt64[]
            end
        else
            println("Error: copy_degree '$copy_degree' not supported! Use either
            'full', 'unit', or 'half'.")
            return
        end
        for smplx_idx in mesh.domains[dom]["simplices"]
            smplx=list_of_simplices[smplx_idx]
            smplx=trace_order[smplx]
            rsmplx=get_reflected_index.(smplx)
            for sector=0:DOS_lim-1
                idx=find_smplx(full_list_of_simplices,get_rotated_index.(smplx, sector))
                ridx=find_smplx(full_list_of_simplices,get_rotated_index.(rsmplx, sector))
                if copy_degree==:full
                    append!(domains[dom]["simplices"],idx)
                    append!(domains[dom]["simplices"],ridx)
                elseif copy_degree==:unit
                    append!(domains[dom*"#$sector"]["simplices"],idx)
                    append!(domains[dom*"#$sector"]["simplices"],ridx)
                elseif copy_degree==:half
                    append!(domains[dom*"#$sector.0"]["simplices"],idx)
                    append!(domains[dom*"#$sector.1"]["simplices"],ridx)
                end
            end
        end
    end
    #create new mesh object
    symmetry_info=(DOS,naxis,nxbloch,nbody,nxsymmetry,n)
    full_mesh=Mesh(mesh.file,fpoints,UInt16[],triangles, tetrahedra, domains, file_txt, zeros(UInt32,length(triangles)), symmetry_info)
    return full_mesh
end
