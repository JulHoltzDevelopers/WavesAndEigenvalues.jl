function readvtu(filename)
    points = []
    tetrahedra = Array{Int64,1}[]
    data=Dict()
    open(filename) do fid
        while !eof(fid)
            line = readline(fid)
            #println(idx,"::",line)
            line=split(line)
            if line[1]=="<Points>"
                line=readline(fid)
                line=split(line,('<','>'))
                points=split(line[3])
                points=parse.(Float64,points)
                points=reshape(points,3,length(points)÷3)
            #     line=split(readline(fid))
            #     while line[1][1]!='<'
            #         numbers=parse.(Float64,line)
            #         append!(points,[numbers[1:3]])
            #         if length(numbers)==6
            #             append!(points,[numbers[4:6]])
            #         end
            #         line=split(readline(fid))
            #     end
            end

            if line[1]=="<Cells>"
                line=readline(fid)
                line=split(line,('<','>'))
                tetrahedra=split(line[3])
                tetrahedra=parse.(UInt32,tetrahedra)
                tetrahedra=reshape(tetrahedra,4,length(tetrahedra)÷4)
                # line=readline(fid)
                # line=split(readline(fid))
                # tet=[]
                # while line[1][1]!='<'
                #     numbers=parse.(Int,line)
                #     while length(numbers)!=0 && length(tet)!=4
                #         append!(tet,popfirst!(numbers))
                #     end
                #     append!(tetrahedra,[tet.+1])
                #     tet=numbers[1:end]
                #     if length(tet)==4
                #         append!(tetrahedra,[tet.+1])
                #         tet=[]
                #     end
                #     line=split(readline(fid))
                # end
            end

            if line[1]=="<PointData>"
                while true
                    line=readline(fid)
                    if strip(line)=="</PointData>"
                        break
                    end
                    line=split(line,('<','>'))
                    name=split(line[2],'"')[4]
                    dates=split(line[3])
                    dates=parse.(Float64,dates)
                    data[name]=dates
                    #read closing DataArray tag
                    line=readline(fid)
                end
            end

        end
    end
    return points, tetrahedra, data
end

points, tetrahedra, data = readvtu("/Users/gmensah/Desktop/shape/eth/eth_lin.vtu")
points, tetrahedra, data_fd = readvtu("/Users/gmensah/Desktop/shape/eth/eth_fd_lin.vtu")

data_diff=Dict()
for key in keys(data)
    data_diff[key]=data[key].-data_fd[key]
end

for val in ("x","y","z")
    idx="discrete adjoint $val[150.34 + 0.0im]Hz"
    maxdat=maximum(abs.(data[idx]))
    maxdiff=maximum(abs.(data_diff[idx]))
    maxrat=abs.(data_diff[idx]./data[idx])
    maxrat=maximum(maxrat[.!isnan.(maxrat)])
    sig_mask=sign.(data[idx]).*sign.(data_fd[idx]).==-1
    sig=any(sig_mask)
    println(any(abs.(data[idx][sig_mask]).>0.02))
    println("max data: $maxdat  max diff:$maxdiff max rat: $maxrat sig: $sig")
end
