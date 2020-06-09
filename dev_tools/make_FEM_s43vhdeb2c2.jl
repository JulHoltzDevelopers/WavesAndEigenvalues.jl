include("make_FEM_header.jl")
## matrices
function local_matrix(order=1,corder=2, crdnt=0, crdnt_adj=0 ;symmetric=false,typ=:none)
    f=F[order]
    if typeof(corder)==SymPy.Sym
        coeff=corder
    elseif corder==0
        coeff=1
    else
        coeff=sum(C[corder])
    end

    if typ==:boundary
        f=f[surf[order]]
    end

    M=Array{Any}(undef,length(f),length(f))
        for (i,fi) in enumerate(f)
            if typ==:nabla
                    fi=grad(fi)
            end
            for (j,fj) in enumerate(f)
                if symmetric && j<i
                    continue
                end
                println("$i $j"); flush(stdout)
                if typ==:nabla
                    fj=grad(fj)
                    intg=matmul(fj,A,fi)
                elseif typ==:boundary
                    intg=(fj.subs(z,0)*fi.subs(z,0)).simplify()
                elseif length(crdnt) == 2
                    intg=fi*SymPy.diff(fj, X[crdnt[1]], X[crdnt[2]])
                elseif  crdnt==0 && crdnt_adj==0
                    intg=(fi*fj).simplify()
                elseif crdnt!=0 && crdnt_adj==0
                    intg=(fi*SymPy.diff(fj,X[crdnt])).simplify()
                elseif crdnt==0 && crdnt_adj!=0
                    intg=(SymPy.diff(fi,X[crdnt_adj])*fj).simplify()
                elseif crdnt!=0 && crdnt_adj!=0
                    intg=(SymPy.diff(fi,X[crdnt_adj])*SymPy.diff(fj,X[crdnt])).simplify()
                end
                if typ==:boundary
                    M[i,j]=tri_integrate((coeff*intg).subs(z,0))
                else
                    M[i,j]=tet_integrate(coeff*intg)
                end
                if symmetric
                    M[j,i]=M[i,j]
                end
            end
        end

    return M
end

##

#Dx=local_matrix(1,diff(sum(C[1])),1)
#Dy=local_matrix(1,0,2)
#Dz=local_matrix(1,0,3)

##
function print_matrix(M,coeff=false,symmetric=false;diff="")
    N=size(M,1)
    if coeff
        txt=""
        for i=1:N
            for j=(symmetric ? i : 1) : N
                txt*="M$diff[$i,$j]="*string(M[i,j])*"\n"
            end
        end
        if symmetric
            for i=1:N
                for j=i+1:N
                    txt*="M$diff[$j,$i]=M$diff[$i,$j]\n"
                end
            end
        end
    else
        txt="["
        for i=1:N
            for j=1:N
                if j!=N
                    txt*=string(M[i,j])*" "
                else
                    txt*=string(M[i,j])
                end
            end
            if i!=N
                txt*=";\n"
            else
                txt*="]"
            end
        end
    end

    txt=replace(txt,"A11"=>"A[1,1]")
    txt=replace(txt,"A12"=>"A[1,2]")
    txt=replace(txt,"A13"=>"A[1,3]")
    txt=replace(txt,"A21"=>"A[2,1]")
    txt=replace(txt,"A22"=>"A[2,2]")
    txt=replace(txt,"A23"=>"A[2,3]")
    txt=replace(txt,"A31"=>"A[3,1]")
    txt=replace(txt,"A32"=>"A[3,2]")
    txt=replace(txt,"A33"=>"A[3,3]")

    txt=replace(txt,"B11"=>"B[1,1]")
    txt=replace(txt,"B12"=>"B[1,2]")
    txt=replace(txt,"B13"=>"B[1,3]")
    txt=replace(txt,"B21"=>"B[2,1]")
    txt=replace(txt,"B22"=>"B[2,2]")
    txt=replace(txt,"B23"=>"B[2,3]")
    txt=replace(txt,"B31"=>"B[3,1]")
    txt=replace(txt,"B32"=>"B[3,2]")
    txt=replace(txt,"B33"=>"B[3,3]")

    for i=1:10
        txt=replace(txt,"c$i^2"=>"cc[$i,$i]")
        txt=replace(txt,"b$i^2"=>"bb[$i,$i]")
        for j=1:10
            txt=replace(txt,"c$i*c$j"=>"cc[$i,$j]")
            txt=replace(txt,"b$i*c$j"=>"bc[$i,$j]")
            txt=replace(txt,"b$i*b$j"=>"bb[$i,$j]")
        end
    end


    println(txt)
    return txt
end
##
#M =local_matrix(1,0) mass matrix
#M=local_matrix(1,0,0,3)
#M=local_matrix(1,diff(sum(C[1]),z),0)
for i=1:3
    for j=1:3
        if j<i
            continue
        end
    M=local_matrix(3,(sum(b[2])*sum(C[2])).expand(),(i,j),typ=:diff,symmetric=false)
    txt=print_matrix(M,true,false,diff="$i$j")
        open("M$i$i.jl","w") do fid
            write(fid,txt)
        end
    end
end
