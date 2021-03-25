export discrete_adjoint_shape_sensitivity, normalize_sensitivity, normal_sensitivity, forward_finite_differences_shape_sensitivity, bound_mass_normalize, blochify_surface_points!
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
    v0Adj/=conj(v0Adj'*L(sol.params[sol.eigval],1)*v0)

    #detect unit mesh
    unit=false
    bloch=false
    if mesh.dos!=1
        unit= mesh.dos.unit
    end
    if unit
        b=:b
    else
        b=:__none__
    end

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
        mesh_h=Mesh("mesh_h",deepcopy(mesh.points),mesh.lines,mesh.triangles,mesh.tetrahedra,domains,"constructed from mesh", mesh.tri2tet,mesh.dos)
        #get cooridnates of current surface point
        pnt=mesh_h.points[:,pnt_idx]

        if unit
            bloch=false
                if 0 < pnt_idx-mesh.dos.naxis <= mesh.dos.nxbloch
                    #identify bloch image point
                    pnt_bloch_idx=size(mesh.points,2)-mesh.dos.nxbloch+(pnt_idx-mesh.dos.naxis)
                    pnt_bloch=mesh.points[:,pnt_bloch_idx]#[:,end-mesh.dos.nxbloch+pnt_idx]
                    bloch=true
                #elseif size(mesh.points,2)-mesh.dos.nxbloch<pnt_idx
                #    #skip analysis if pnt is bloch image point
                #    continue
                elseif pnt_idx<=mesh.dos.naxis
                    #skip analysis if pnt is axis pnt
                    if output
                        ProgressMeter.next!(p)
                    end
                    continue
                end
        end

        #perturb coordinates for all three directions and compute operator derivative by central FD
        #then use operator derivative to compute eigval sensitivity by adjoint approach
        for crdnt=1:3
            mesh_h.points[:,pnt_idx]=pnt

            if unit
                X=get_cylindrics(pnt)
                mesh_h.points[:,pnt_idx].+=h.*X[:,crdnt]
                if bloch
                   mesh_h.points[:,pnt_bloch_idx]=pnt_bloch
                   X_bloch=get_cylindrics(pnt_bloch)
                   #println("####X: $pnt_idx")
                   mesh_h.points[:,pnt_bloch_idx].+=h.*X_bloch[:,crdnt]
                end
            else
                mesh_h.points[crdnt,pnt_idx]+=h
            end
            D_right=discretize(mesh_h, dscrp, C, mass_weighting=false,b=b)
            if unit
                mesh_h.points[:,pnt_idx].-=2*h.*X[:,crdnt]
                if bloch
                   mesh_h.points[:,pnt_bloch_idx].-=2*h.*X_bloch[:,crdnt]
                end
            else
                mesh_h.points[crdnt,pnt_idx]-=2h
            end
            D_left=discretize(mesh_h, dscrp, C, mass_weighting=false,b=b)
            if unit
                D_right.params[b]=1
                D_left.params[b]=1
            end
            D=(D_right(ω0)-D_left(ω0))/(2*h)
            #println("#######")
            #println("size::$(size(D))")
            #println("###D:$(-v0Adj'*D*v0)   idx:$pnt_idx")
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
    nsens=bound_mass_normalize(surface_points,normal_vectors,tri_mask,sens)

Normalize shape sensitivity with boundary mass matrix.
"""
function bound_mass_normalize(surface_points,normal_vectors,tri_mask,mesh,sens)
    M=[1/12 1/24 1/24;
       1/24 1/12 1/24;
       1/24 1/24 1/12]

       II=[]
       JJ=[]
       MM=[]

       for (idx,tri) in enumerate(mesh.triangles)
           ii=[tri[1] tri[2] tri[3];
           tri[1] tri[2] tri[3];
           tri[1] tri[2] tri[3]]
           jj=ii'
           mm=M.*LinearAlgebra.norm(normal_vectors[:,idx])
           append!(II,ii[:])
           append!(JJ,jj[:])
           append!(MM,mm[:])
       end

       #relabel points
       D=Dict()
       for (idx,pnt) in enumerate(surface_points)
           D[pnt]=idx
       end
       for jdx in 1:length(II)
            II[jdx]=D[II[jdx]]
            JJ[jdx]=D[JJ[jdx]]
       end
    B=SparseArrays.sparse(II,JJ,MM)
    B=SparseArrays.lu(B)
    nsens=zeros(ComplexF64,size(sens))
    for i=1:3
        nsens[i,surface_points]=B\sens[i,surface_points]
    end
    return nsens
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
    #check whether is unit mesh
    unit=0
    if mesh.dos!=1
        unit=mesh.dos.unit
    end

    sens=zeros(ComplexF64,size(mesh.points))

    if unit!=0
        n_iterations=length(surface_points)-mesh.dos.nxbloch
    else
        n_iterations=length(surface_points)
    end

    if output
        p = ProgressMeter.Progress(n_iterations,desc="FD Sensitivity... ", dt=1,
             barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=20)
    end



    for idx=1:n_iterations
        #global G
        pnt_idx=surface_points[idx]
        domains=assemble_connected_domain(idx,mesh::Mesh,dscrp,tri_mask,tet_mask)


        #get coordinates of current surface point
        pnt=mesh.points[:,pnt_idx]
        for crdnt=1:3
            #construct mesh with reduced domains
            mesh_h=Mesh("mesh_h",deepcopy(mesh.points),mesh.lines,mesh.triangles,mesh.tetrahedra,domains,"constructed from mesh", mesh.tri2tet,mesh.dos)
            mesh_hl=deepcopy(mesh_h)
            #forward finite difference
            #mesh_h.points[crdnt,pnt_idx]+=h
            #G=discretize(mesh_h, dscrp, C)
            #D_center=discretize(mesh_h, dscrp, C, mass_weighting=true)
            if unit>0
                X=get_cylindrics(pnt)
                mesh_h.points[:,pnt_idx].+=h.*X[:,crdnt]
                mesh_hl.points[:,pnt_idx].-=h.*X[:,crdnt]
                #bidx=pnt_idx-mesh.dos.naxis #(potential) bloch point index
                #if 0 < bidx <= mesh.dos.nxbloch
                if 0 < pnt_idx-mesh.dos.naxis <= mesh.dos.nxbloch
                    #identify bloch image point
                    bidx=size(mesh.points,2)-mesh.dos.nxbloch+(pnt_idx-mesh.dos.naxis)
                    bpnt=mesh.points[:,bidx]
                    bX=get_cylindrics(bpnt)
                    mesh_h.points[:,bidx].+=h.*bX[:,crdnt]
                    mesh_hl.points[:,bidx].-=h.*bX[:,crdnt]
                    # wrong
                    # mesh_h.points[:,bidx]=mesh_h.points[:,pnt_idx]
                    # mesh_h.points[2,bidx]*=-1 #wrong
                else
                    bidx=0 #sentinel value
                end
            else
                mesh_h.points[crdnt,pnt_idx]+=h
                mesh_hl.points[crdnt,pnt_idx]-=h
            end
            if unit>0
                D_right=discretize(mesh_h, dscrp, C, mass_weighting=true,b=:b)
                D_left=discretize(mesh_hl, dscrp, C, mass_weighting=true,b=:b)
            else
                D_right=discretize(mesh_h, dscrp, C, mass_weighting=true)
                D_left=discretize(mesh_hl, dscrp, C, mass_weighting=true)
            end
            #G=deepcopy(L)
            G=LinearOperatorFamily(["ω","λ"],complex([0.,Inf]))
            for (key,val) in L.params
                G.params[key]=val
            end
            for idx in 1:length(D_left.terms)
                symbol=L.terms[idx].symbol
                if symbol!="__aux__"
                    coeff=L.terms[idx].coeff+D_right.terms[idx].coeff-D_left.terms[idx].coeff
                else
                    coeff=L.terms[idx].coeff
                end
                func=L.terms[idx].func
                operator=L.terms[idx].operator
                params=L.terms[idx].params

                push!(G,Term(coeff,func,params,symbol,operator))
            end
            new_sol, n, flag = householder(G,sol.params[sol.eigval],maxiter=5, output = false, nev=3,order=3)
            sens[crdnt,pnt_idx]=(new_sol.params[new_sol.eigval]-sol.params[sol.eigval])/(2h)
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
##
function blochify_surface_points!(mesh::Mesh, surface_points, tri_mask, tet_mask)
    for idx in surface_points
        bidx=idx -mesh.dos.naxis
        if 0 < bidx <= mesh.dos.nxbloch
            append!(tri_mask[idx],tri_mask[end-mesh.dos.nxbloch+bidx])
            append!(tet_mask[idx],tet_mask[end-mesh.dos.nxbloch+bidx])
            unique!(tri_mask)
            unique!(tet_mask)
        end
    end
    return nothing
end

##
"helper function to get local cylindric basis vectors"
function get_cylindrics(pnt)
    #r,phi,z=X
    X=zeros(3,3)
    X[:,3]=[0;0;1]
    X[:,1]=copy(pnt)
    X[3,1]=0.0
    X[:,1]./=LinearAlgebra.norm(X[:,1])
    X[:,2]=LinearAlgebra.cross(X[:,3],X[:,1])
    return X
end
