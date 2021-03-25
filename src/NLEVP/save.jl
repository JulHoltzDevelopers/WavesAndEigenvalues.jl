import Dates
function save(fname,L::Solution)
  date=string(Dates.now(Dates.UTC))
  open(fname,"w") do f
    write(f,"# Solution version 0\n")
    write(f,"#"*date*"\n")
    write(f,"params=[")
    for (key,value) in L.params
      write(f,"(:$(string(key)),$value),\n")
    end
    write(f,"]\n")
    write(f,"eigval=:$(L.eigval)\n")

    write(f,"v=")
    vector_write(f,L.v)
    write(f,"]\n")

    write(f,"v_adj=")
    vector_write(f,L.v_adj)

    write(f,"[eigval_pert]\n")
    for (key,value) in L.eigval_pert
      write(f,"\t[eigval_pert.$key]\n")
      if typeof(value)<:Tuple
        write(f,"\t\tnum=")
        vector_write(f,value[1])
        write(f,"\t\tden=")
        vector_write(f,value[2])
      else
        write(f,"\t\tnum=")
        vector_write(f,value)
      end
    end

    write(f,"[v_pert]\n")
    for (key,value) in L.v_pert
      write(f,"\t[v_pert.$key]\n")
      if typeof(value)<:Tuple
        write(f,"\t\t[v_pert.$key.num]\n")
        for (idx,val) in enumerate(value[1])
          write(f,"\t\t\t[v_pert.$key.num.$idx]\n")
          write(f,"\t\t\tv=")
          vector_write(f,val)
        end

        write(f,"\t\t[v_pert.$key.den]\n")
        for (idx,val) in enumerate(value[2])
          write(f,"\t\t\t[v_pert.$key.den.$idx]\n")
          write(f,"\t\t\tv=")
          vector_write(f,val)
        end
      else
        write(f,"\t\t[v_pert.$key.num]\n")
        for (idx,val) in enumerate(value)
          write(f,"\t\t\t[v_pert.$key.num.$idx]\n")
          write(f,"\t\t\tv=")
          vector_write(f,val)
        end
      end
    end
    # for (idx,val) in enumerate(L.v_pert)
    #   write(f,"\t[v_pert.$idx]\n")
    #   write(f,"\t\tv=")
    #   vector_write(f,val)
    # end
  end
end
#TODO: save v_pert
function vector_write(f,V)
  write(f,"[")
  for v in V #TODO: mehr Spalten
    if imag(v)>=0
      write(f,"$(real(v))+$(imag(v))im,")
    else
      write(f,"$(real(v))$(imag(v))im,")
    end
  end
  write(f,"]\n")
  return
end



#include("toml.jl")

function read_sol(f)
    D=read_toml(f)
    eigval_pert=Dict{Symbol,Any}()
    for (key,value) in D["/eigval_pert"]
      if haskey(value,"den")
        eigval_pert[Symbol(key[2:end])]=(value["num"],value["den"])
      else
        eigval_pert[Symbol(key[2:end])]=value["num"]
      end
    end
    params=Dict{Symbol,Complex{Float64}}()
    for (sym,val) in D["params"]
      params[sym]=val
    end
    eigval=D["eigval"]
    v=D["v"]
    v_adj=D["v_adj"]


    v_pert=Dict{Symbol,Any}()
    for (key,value) in D["/v_pert"]
      if haskey(value,"/den")
        N = length(value["/den"])
        B=Array{Array{Complex{Float64}},1}(undef,N)
        for idx=1:N
          B[idx] = value["/den"]["/$idx"]["v"]
        end
        N = length(value["/num"])
        A=Array{Array{Complex{Float64}},1}(undef,N)
        for idx=1:N
          A[idx] = value["/num"]["/$idx"]["v"]
        end
        v_pert[Symbol(key[2:end])]=A,B
      else
        N=length(value["/num"])
        V=Array{Array{Complex{Float64}},1}(undef,N)
        for idx=1:N
          V[idx] = value["/num"]["/$idx"]["v"]
        end
        v_pert[Symbol(key[2:end])]=V
      end
    end
    #TODO: proper type initialization
    # N=length(D["/v_pert"])
    # v_pert=Array{Array{Complex{Float64}},1}(undef,N)
    # for idx=1:N
    #   v_pert[idx]=D["/v_pert"]["/$idx"]["v"]
    # end
    return Solution(params,v,v_adj,eigval,eigval_pert,v_pert)
end
