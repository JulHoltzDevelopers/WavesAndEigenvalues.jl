import WavesAndEigenvalues.Meshutils: insert_smplx!, find_smplx,sort_smplx
import SparseArrays, LinearAlgebra
##

##
B = bound_mass(surface_points,normal_vectors,tri_mask)
B=SparseArrays.lu(B)
##
nsens=zeros(ComplexF64,size(sens))
for i=1:3
    nsens[i,surface_points]=B\sens[i,surface_points]
end
##
function mass_normalize_sensitivity(surface_points,normal_vectors,tri_mask,sens)
    #preallocate arrays
    A=Array{Float64}(undef,length(normal_vectors))#
    V=Array{Float64}(undef,length(normal_vectors))
    normed_sens=zeros(ComplexF64,size(normal_vectors))
    for (crdnt,v) in enumerate(([1.0; 0.0; 0.0],[0.0; 1.0; 0.0],[0.0; 0.0; 1.0]))
        #compute triangle area
        for idx =1:size(normal_vectors,2)
            A[idx]=LinearAlgebra.norm(normal_vectors[:,idx])/2
            #V[idx]=LinearAlgebra.norm(LinearAlgebra.cross(normal_vectors[:,idx],v))/6
            #V[idx]=abs(LinearAlgebra.dot(normal_vectors[:,idx],v))/6
        end

        for idx=1:length(surface_points)
            pnt=surface_points[idx]
            tris=tri_mask[idx]
            vol=sum(abs.(V[tris]))#common volume flow of all adjacent triangles
            vol=0.0
            for tri in tris
                vol+=LinearAlgebra.dot(normal_vectors[:,idx],v)/6
            end
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
