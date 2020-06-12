include("make_FEM_header.jl")
###
function local_vector(order=1,corder=2, crdnt=0, crdnt_adj=0, typ=:none )
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

    if typ==:gradient
        M=Array{Any}(undef,length(f),3)
    else
        M=Array{Any}(undef,length(f),1)
    end
    for (i,fi) in enumerate(f)
        if typ==:boundary
            M[i,1]=tri_integrate((fi*coeff).subs(z, 0).simplify()).simplify()
        elseif typ==:gradient
            M[i,:]=grad(fi*coeff)
        else
            M[i,1]=tet_integrate(fi*coeff).simplify()
        end
    end

    return M
end

##
function print_matrix(M,coeff=false,symmetric=false;diff="")
    N=size(M)
    if coeff
        txt=""
        for i=1:N[1]
            for j=(symmetric ? i : 1) : N[2]
                txt*="M$diff[$i,$j]="*string(M[i,j])*"\n"
            end
        end
        if symmetric
            for i=1:N[1]
                for j=i+1:N[2]
                    txt*="M$diff[$j,$i]=M$diff[$i,$j]\n"
                end
            end
        end
    else
        txt="["
        for i=1:N[1]
            for j=1:N[2]
                if j!=N[2]
                    txt*=string(M[i,j])*" "
                else
                    txt*=string(M[i,j])
                end
            end
            if i!=N[1]
                txt*=";\n"
            else
                txt*="]"
            end
        end
    end



    for i=1:3
        for j=1:3
            txt=replace(txt,"A$i$j"=>"A[$i,$j]")
            txt=replace(txt,"B$i$j"=>"B[$i,$j]")
        end
    end


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




M=local_vector(3,0,0,0,:vol)

txt=print_matrix(M,false)

txt=replace(txt,";"*'\n'=>" ")
println(txt)
