mesh=Mesh("./examples/tutorials/Rijke_mm.msh",scale=0.001)

function helmhol(mesh,order)
    triangles,tetrahedra=aggregate_elements(mesh,order)
    function stiff(J)
        if order==:1
            return -s43nv1nu1(J)
        elseif order==:2
            return -s43nv2nu2(J)
        elseif order==:h
            return -s43nvhnuh(J)#stiffh(J)
        end
    end
    function mass(J)
        if order==:1
            return s43v1u1(J)
        elseif order==:2
            return s43v2u2(J)
        elseif order==:h
            return s43vhuh(J)#massh(J)
        end
    end

    L=LinearOperatorFamily(["ω","λ"],complex([0.,Inf]))

    MM=ComplexF64[]
    II=Int64[]
    JJ=Int64[]
    for smplx in tetrahedra
            J=CooTrafo(mesh.points[:,smplx[1:4]])
            mm=347.0^2*stiff(J)
            ii,jj=create_indices(smplx)
            append!(MM,mm[:])
            append!(II,ii[:])
            append!(JJ,jj[:])
    end

    M=SparseArrays.sparse(II,JJ,MM)
    push!(L,Term(M,(),(),"","K"))

    MM=ComplexF64[]
    II=Int64[]
    JJ=Int64[]
    for smplx in tetrahedra
            J=CooTrafo(mesh.points[:,smplx[1:4]])
            mm=mass(J)
            ii,jj=create_indices(smplx)
            append!(MM,mm[:])
            append!(II,ii[:])
            append!(JJ,jj[:])
    end

    M=SparseArrays.sparse(II,JJ,MM)
    push!(L,Term(M,(pow2,),((:ω,),),"ω^2","M"))
    push!(L,Term(-M,(pow1,),((:λ,),),"-λ","__aux__"))
    return L
end
##
L=helmhol(mesh,:h)
##
sol,nn,flag=householder(L,340*2*pi,maxiter=5,output=true)
##
c=ones(length(mesh.tetrahedra))*347
c=ones(size(mesh.points,2))*347

dscrp=Dict()
dscrp["Interior"]=(:interior,())
H=WavesAndEigenvalues.Helmholtz.discretize(mesh,dscrp,c,el_type=2,c_type=1)
sol,nn,flag=householder(H,330*2*pi,maxiter=2,output=true)
###

X=[0 1 0 0;
   1 0 0 1;
   0 0 1 1.]
X=rand(3,4)*100
J=CooTrafo(X)
m1=s43vhuh(J)
m2=massh(J)
##
