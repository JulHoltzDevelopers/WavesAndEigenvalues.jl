
function readvtu(filename)
    points = []
    tetrahedra = Array{Int64,1}[]

    open(filename) do fid
        while !eof(fid)
            line = readline(fid)
            #println(idx,"::",line)
            line=split(line)
            if line[1]=="<Points>"
                line=readline(fid)
                line=split(readline(fid))

                while line[1][1]!='<'
                    numbers=parse.(Float64,line)
                    append!(points,[numbers[1:3]])
                    if length(numbers)==6
                        append!(points,[numbers[4:6]])
                    end
                    line=split(readline(fid))
                end
            end

            if line[1]=="<Cells>"
                line=readline(fid)
                line=split(readline(fid))
                tet=[]
                while line[1][1]!='<'
                    numbers=parse.(Int,line)
                    while length(numbers)!=0 && length(tet)!=4
                        append!(tet,popfirst!(numbers))
                    end
                    append!(tetrahedra,[tet.+1])
                    tet=numbers[1:end]
                    if length(tet)==4
                        append!(tetrahedra,[tet.+1])
                        tet=[]
                    end
                    line=split(readline(fid))
                end
            end

        end
    end
    return points, tetrahedra
end

points, tetrahedra = readvtu("mesh.vtu")

import ProgressMeter

p = ProgressMeter.Progress(length(tetrahedra),desc="Assemble all triangles ", dt=1,
     barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
     barlen=20)
triangles=[]
simplices=[]
for (idx,tet) in enumerate(tetrahedra)
    append!(simplices,idx)

    for tri in (tet[[1,2,3]],tet[[1,2,4]],tet[[1,3,4]],tet[[2,3,4]])
        idx,ins=sort_smplx(triangles,tri)
        if ins
            insert!(triangles,idx,tri)
        else
            deleteat!(triangles,idx)
        end
    end
    ProgressMeter.next!(p)
end
# convert points

Points=zeros(Float64,3,length(points))
for (idx,pnt) in enumerate(points)
    Points[:,idx]=pnt
end
#
x_min=minimum(Points[1,:])
x_max=maximum(Points[1,:])
#
domains=Dict()
domain=Dict()
domain["simplices"]=simplices
domain["dimension"]=0x00000003
domains["Interior"]=domain
#
x_min=minimum(Points[1,:])
x_max=maximum(Points[1,:])
simplices=[]
for (tri_idx,tri) in enumerate(triangles)
    if all(abs.(Points[1,tri].-x_max).<=0.0001)
        insert_smplx!(simplices,tri_idx)
    end
end
domain=Dict()
domain["simplices"]=simplices
domain["dimension"]=0x00000002
domains["Outlet"]=domain
#
tri2tet=zeros(UInt32,length(triangles))
tri2tet[1]=0xffffffff

mesh=Mesh("mesh.vtu",Points,[],triangles,tetrahedra,domains,"mesh.vtu",tri2tet,1)
vtk_write("eth_diff", mesh, data_diff)
