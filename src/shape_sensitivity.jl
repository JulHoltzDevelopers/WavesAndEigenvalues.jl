export discrete_adjoint_shape_sensitivity, normalize_sensitivity, normal_sensitivity, forward_finite_differences_shape_sensitivity
"""
    sens=discrete_adjoint_shape_sensitivity(mesh::Mesh,dscrp,C,surface_points,tri_mask,tet_mask,L,sol;h=1E-9,output=true)

Compute shape sensitivity of mesh `mesh`

#Arguments
- `mesh::Mesh`: Mesh
-  ...

#Notes
The method is based on a discrete adjoint approach. The derivative of the
operator with respect to a point displacement is, however, computed by central
finite differences.
"""
function discrete_adjoint_shape_sensitivity(mesh::Mesh,dscrp,C,surface_points,tri_mask,tet_mask,L,sol;h=1E-9,output=true)
    ω0=sol.params[sol.eigval]

    ## normalize base
    #TODO: make this normalization standard to avoid passing L
    v0=sol.v
    v0Adj=sol.v_adj
    v0/=sqrt(v0'*v0)
    v0Adj/=v0Adj'*L(sol.params[sol.eigval],1)*v0



    sens=zeros(ComplexF64,3,size(mesh.points,2))
    if output
        p = ProgressMeter.Progress(length(surface_points),desc="DA Sensitivity... ", dt=1,
             barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=20)
    end
    for idx=1:length(surface_points)
        #global D, D_left, D_right, domains
        pnt_idx=surface_points[idx]
        domains=Dict()
        #shorten domain simplex-list s.t. only simplices that are connected to current
        #surface point occur.
        for dom in keys(dscrp)
            domain=deepcopy(mesh.domains[dom])
            dim=domain["dimension"]
            simplices=[]
            if dim==2
                smplx_list=tri_mask[idx]
            elseif dim==3
                smplx_list=tet_mask[idx]
            else
                smplx_list=[]
            end

            for smplx_idx in domain["simplices"]
                if smplx_idx in smplx_list
                    append!(simplices,smplx_idx)
                end
            end
            domain["simplices"]=simplices
            domains[dom]=domain
        end
        #construct mesh with reduced domains
        mesh_h=Mesh("mesh_h",deepcopy(mesh.points),mesh.lines,mesh.triangles,mesh.tetrahedra,domains,"constructed from mesh", mesh.tri2tet,[])
        #get cooridnates of current surface point
        pnt=mesh_h.points[:,pnt_idx]
        #perturb coordinates for all three directions and compute operator derivative by central FD
        #then use operator derivative to compute eigval sensitivity by adjoint approach
        for crdnt=1:3
            mesh_h.points[:,pnt_idx]=pnt
            mesh_h.points[crdnt,pnt_idx]+=h
            D_right=discretize(mesh_h, dscrp, C, mass_weighting=false)
            mesh_h.points[crdnt,pnt_idx]-=2h
            D_left=discretize(mesh_h, dscrp, C, mass_weighting=false)
            D=(D_right(ω0)-D_left(ω0))/(2*h)
            sens[crdnt,pnt_idx]=-v0Adj'*D*v0
        end

        if output
            ProgressMeter.next!(p)
        end
    end
    return sens
end


"""
    normed_sens=normalize_sensitivity(surface_points,normal_vectors,tri_mask,sens)

Normalize shape sensitivity with directed area of adjacent triangles.
"""
function normalize_sensitivity(surface_points,normal_vectors,tri_mask,sens)
    #preallocate arrays
    A=Array{Float64}(undef,length(normal_vectors))#
    V=Array{Float64}(undef,length(normal_vectors))
    normed_sens=zeros(ComplexF64,size(normal_vectors))
    for (crdnt,v) in enumerate(([1.0; 0.0; 0.0],[0.0; 1.0; 0.0],[0.0; 0.0; 1.0]))
        #compute triangle area and volume flow
        for idx =1:size(normal_vectors,2)
            A[idx]=LinearAlgebra.norm(normal_vectors[:,idx])/2
            #V[idx]=LinearAlgebra.norm(LinearAlgebra.cross(normal_vectors[:,idx],v))/6
            V[idx]=abs(LinearAlgebra.dot(normal_vectors[:,idx],v))/6
        end

        for idx=1:length(surface_points)
            pnt=surface_points[idx]
            tris=tri_mask[idx]
            vol=sum(abs.(V[tris]))#common volume flow of all adjacent triangles
            if vol==0 #TODO: compare against some reference value to check for numerical zero
                #println("$idx")
                continue
            end
            for tri in tris
                weight=abs(V[tri])/vol
                if A[tri]>0
                    normed_sens[crdnt,tri]+=sens[crdnt,pnt]/A[tri]*weight
                end
                if isnan(normed_sens[crdnt,idx])
                    println("$idx: $(vol==0)")
                    println("idx:$idx weight:$weight A:$(A[tri]) vol:$(vol==0)")
                end
            end
        end
    end
    return normed_sens
end

"""
    normal_sensitivity(normal_vectors, normed_sens)

Convert surface normalized shape_gradient `normed_sens` to surface_normal shape
sensitivity.
"""
function normal_sensitivity(normal_vectors, normed_sens)
    normal_sens=Array{ComplexF64}(undef,size(normal_vectors,2))
    for idx =1:size(normal_vectors,2)
        n=normal_vectors[:,idx]
        s=normed_sens[:,idx]
        v=n./LinearAlgebra.norm(n)
        normal_sens[idx]=LinearAlgebra.dot(v,s)
    end
    return normal_sens
end

## FD approach
function forward_finite_differences_shape_sensitivity(mesh::Mesh,dscrp,C,surface_points,tri_mask,tet_mask,L,sol;h=1E-9,output=true)
    sens=zeros(ComplexF64,size(mesh.points))
    if output
        p = ProgressMeter.Progress(length(surface_points),desc="FD Sensitivity... ", dt=1,
             barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=20)
    end
    for idx=1:length(surface_points)
        #global G
        pnt_idx=surface_points[idx]
        domains=assemble_connected_domain(idx,mesh::Mesh,dscrp,tri_mask,tet_mask)

        #construct mesh with reduced domains
        mesh_h=Mesh("mesh_h",deepcopy(mesh.points),mesh.lines,mesh.triangles,mesh.tetrahedra,domains,"constructed from mesh", mesh.tri2tet,[])
        #get cooridnates of current surface point
        pnt=mesh_h.points[:,pnt_idx]

        for crdnt=1:3
        #forward finite difference
            #mesh_h.points[crdnt,pnt_idx]+=h
            #G=discretize(mesh_h, dscrp, C)
            D_center=discretize(mesh_h, dscrp, C, mass_weighting=true)
            mesh_h.points[crdnt,pnt_idx]+=h
            D_right=discretize(mesh_h, dscrp, C, mass_weighting=true)
            G=deepcopy(L)
            G=LinearOperatorFamily(["ω","λ"],complex([0.,Inf]))
            for (key,val) in L.params
                G.params[key]=val
            end
            for idx in 1:length(D_center.terms)
                coeff=L.terms[idx].coeff+D_right.terms[idx].coeff-D_center.terms[idx].coeff
                func=L.terms[idx].func
                operator=L.terms[idx].operator
                params=L.terms[idx].params
                symbol=L.terms[idx].symbol
                push!(G,Term(coeff,func,params,symbol,operator))
            end
            new_sol, n, flag = householder(G,sol.params[sol.eigval],maxiter=4, output = false, n_eig_val=3)
            sens[crdnt,pnt_idx]=(new_sol.params[new_sol.eigval]-sol.params[sol.eigval])/(h)
        end

        if output
            #update progressMeter
            ProgressMeter.next!(p)
        end
    end
    return sens
end
##
function assemble_connected_domain(idx,mesh::Mesh,dscrp,tri_mask,tet_mask)
    domains=Dict()
    for dom in keys(dscrp)
        domain=deepcopy(mesh.domains[dom])
        dim=domain["dimension"]
        simplices=[]
        if dim==2
            smplx_list=tri_mask[idx]
        elseif dim==3
            smplx_list=tet_mask[idx]
        else
            smplx_list=[]
        end

        for smplx_idx in domain["simplices"]
            if smplx_idx in smplx_list
                append!(simplices,smplx_idx)
            end
        end
        domain["simplices"]=simplices
        domains[dom]=domain
    end
    return domains
end
