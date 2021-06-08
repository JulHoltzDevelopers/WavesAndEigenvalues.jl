"""
    p=get_p(mesh,v,point,tet_idx=0)

Get pressure value `p` at point `point`in solution vector `v`.
The optional argument tet_idx accelarates finding the tetrahedron enclosing the point.
"""
function get_p(mesh::Mesh,v,point,tet_idx::Int=0)
    CT, order, smplx = initialize_CT(mesh,v,point,tet_idx)

    if order==:lin
        f=f1

    elseif order==:quad
        f=f2

    elseif order==:herm
        f=fh
    else
        return nothing
    end

    return sum(f(CT,point).*v[smplx])
end
"""
    get_n_grad_p(mesh::Mesh,v,point,n,tet_idx::Int=0)

Get directional derivative `dp` of teh pressure w.r.t. the direction `n` at
point `point`in solution vector `v`. The optional argument tet_idx accelarates
finding the tetrahedron enclosing the point.
"""
function get_n_grad_p(mesh::Mesh,v,point,n,tet_idx::Int=0)
    CT,order,smplx = initialize_CT(mesh,v,point,tet_idx)

    function gradsrc(J,n,x)
        if order==:lin
            return s43nv1rx(J,n,x)
        elseif order==:quad
            return s43nv2rx(J,n,x)
        elseif order==:herm
            return s43nvhrx(J,n,x)
        end
    end

    return sum(gradsrc(CT,n,point).*v[smplx])
end

function initialize_CT(mesh,v,point,tet_idx=0)
    if tet_idx==0
        tet_idx=find_tetrahedron_containing_point(mesh,point)
    end
    dim=length(v) #TODO:Blochify
    N_points=size(mesh.points,2)
    N_lines=length(mesh.lines)
    N_triangles=length(mesh.triangles)+length(mesh.int_triangles)
    if dim == N_points
        el_type=:lin
    elseif dim == N_points+N_lines
        el_type=:quad
    elseif dim == 4*N_points+N_triangles
        el_type=:herm
    else
        el_type=:unknown
    end
    smplx=aggregate_element(mesh,tet_idx, el_type)
    CT=CooTrafo(mesh.points[:,smplx[1:4]])

    return CT, el_type, smplx
end

##

function aggregate_element(mesh::Mesh, tet_idx, el_type=:lin)
    N_points=size(mesh.points)[2]
    if (el_type in (:quad,:herm) ) &&  length(mesh.lines)==0
        collect_lines!(mesh)
    end
    smplx=mesh.tetrahedra[tet_idx]
    if el_type==:lin
        tet=copy(smplx)

    elseif el_type==:quad
        tet=Array{UInt32}(undef,10)
        tet[1:4]=copy(smplx[:])
        tet[5]=get_line_idx(mesh,smplx[[1,2]])+N_points#find_smplx(mesh.lines,smplx[[1,2]])+N_points #TODO: type stability
        tet[6]=get_line_idx(mesh,smplx[[1,3]])+N_points
        tet[7]=get_line_idx(mesh,smplx[[1,4]])+N_points
        tet[8]=get_line_idx(mesh,smplx[[2,3]])+N_points
        tet[9]=get_line_idx(mesh,smplx[[2,4]])+N_points
        tet[10]=get_line_idx(mesh,smplx[[3,4]])+N_points

    elseif el_type==:herm
        tet=Array{UInt32}(undef,20)
        tet[1:4] = smplx[:].+0
        tet[5:8] = smplx[:].+N_points
        tet[9:12] = smplx[:].+2*N_points
        tet[13:16]= smplx[:].+3*N_points
        for (jdx,tria) in enumerate([smplx[[2,3,4]],smplx[[1,3,4]],smplx[[1,2,4]],smplx[[1,2,3]]])
            fcidx     =  find_smplx(mesh.triangles,tria)
            if fcidx !=0
                tet[16+jdx]   =  fcidx+4*N_points
            else
                fcidx=find_smplx(mesh.int_triangles,tria)
                if fcidx==0
                    println("Error, face not found!!!")
                    return nothing
                end
                tet[16+jdx]   =  fcidx+4*N_points+length(mesh.triangles)
            end
        end

    else
        println("Error: element order $(:el_type) not defined!")
        return nothing
    end

    return tet

end
