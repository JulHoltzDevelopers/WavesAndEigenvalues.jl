"""
Module containing functionality to read and process tetrahedral meshes from gmsh.
"""
module Meshutils
import LinearAlgebra
include("Meshutils_exports.jl")

"""
Definition of the mesh type

# Fields
- `name::String`: the name of the mesh.
- `points::Array`: 3×N array containing the coordinates of the N points defining the mesh.
- `lines::List`: List of simplices defining the edges of the mesh
- `triangles::List`: List of simplices defining the surface triangles of the mesh
- `tetrahedra::List`: List of simplices defining the tetrahedra of the mesh
- `domains::Dict`: Dictionairy defining the domains of the mesh. See comments below.
- `file::String`: path to the file containing the mesh.
- `tri2tet::Array`: Array of length `length(tetrahedra)` containing the indices of the connected tetrahedra.
- `dos`: special field meant to contain symmetry information of highly symmetric meshes.

# Notes
The meshes are supposed to be tetrahedral. All simplices (lines, triangles, and tetrahedra) are stored as list of simplices.
Simplices are lists of integers containing the indices of the points (i.e. the column number in the `points` array) forming the simplex.
This means a line is a two-entry list, a triangle a three-entry list, and a tetrahedron a four entry lists.
For convenience certain entities of the mesh can be further defined in the `domains`dictionary. Each key defines a domain and maps to another dictionary.
This second-level dictionary contains at least two keys: `"dimension"` mapping to the dimension of the specified domain (1,2, or 3) and
`"simplices"` containing a list of integers mapping into the respective simplex lists. More keys may be added to the dictionary to
define additional and/or custom information on the domain. For instance the `compute_volume!` function adds an entry with the domain size.
"""
struct Mesh
    name
    points
    lines
    triangles
    tetrahedra
    domains
    file
    tri2tet
    dos
end
## constructor
"""
    mesh=Mesh(file_name::String; scale=1)

read a tetrahedral mesh from gmsh file into `mesh`.
The optional scaling factor `scale` may be used to scale the units of the mesh.
"""
function Mesh(file_name::String;scale=1)
    #mesh constructor
    #check file format
    local line
    open(file_name,"r") do fid
        line=readline(fid)
        line=readline(fid)
    end
    if line[1]=='4'
        points, lines, triangles, tetrahedra, domains =read_msh4(file_name)
    else
        println("Please, convert msh-file to gmsh's .msh-format version 4!")
        #points, lines, triangles, tetrahedra, domains =read_msh2(file_name)
    end

    #sometimes elements are multiply defined, make them unique
    #TODO: This can be done more memory efficient by reading the msh-file twice
    ulines=Array{UInt32,1}[]
    utriangles=Array{UInt32,1}[]
    utetrahedra=Array{UInt32,1}[]

    for ln in lines
        insert_smplx!(ulines,ln)
    end
    for tri in triangles
        insert_smplx!(utriangles,tri)
    end
    for tet in tetrahedra
        insert_smplx!(utetrahedra,tet)
    end

    for dom in keys(domains)
        for (idx,smplx_idx) in enumerate(domains[dom]["simplices"])
            if domains[dom]["dimension"]==1
                new_idx=find_smplx(ulines,lines[smplx_idx])
            elseif domains[dom]["dimension"]==2
                new_idx=find_smplx(utriangles,triangles[smplx_idx])
            elseif domains[dom]["dimension"]==3
                new_idx=find_smplx(utetrahedra,tetrahedra[smplx_idx])
            end
            domains[dom]["simplices"][idx]=new_idx
        end
        unique!(domains[dom]["simplices"])
    end


    tri2tet=zeros(UInt32,length(utriangles))
    Mesh(file_name, points*scale, ulines, utriangles, utetrahedra, domains, file_name,tri2tet,1)
end

#show mesh
import Base.string
function string(mesh::Mesh)
    txt="""mesh: $(mesh.file)
    #################
    points:     $(size(mesh.points)[2])
    lines:      $(length(mesh.lines))
    triangles:  $(length(mesh.triangles))
    tetrahedra: $(length(mesh.tetrahedra))
    #################
    domains: """
    for key in sort([key for key in keys(mesh.domains)])
        txt*=string(key)*", "
    end
    return txt[1:end-2]
end
import Base.show
function show(io::IO,mesh::Mesh)
  print(io,string(mesh))
end

##
# import Base.getindex(mesh::Mesh, i::String)
#     #parse for sectors
#     splt=split(i,"_#")
#     sctr,hlf=-1,-1
#     dom=splt[1] #domain name
#     #TODO: error message of split has more than 3 entries
#     if length(splt)>1
#         splt=split(splt[2],".")
#         sctr=splt[1]
#         if length(splt)>1 #TODO: length check, see above
#             hlf=splt[2] #TODO: check that this value should be 0 or 1 (assert macro?)
#         end
#     end
#
#     #just call domain
#     domain = Dict()
#     domain["dimension"]=mesh.domains[dom]["dimension"]
#     smplx_list=[]
#     #TODO: compute volume if existent
#     if hlf>-1 #build half-cell
#         for smplx in mesh.domains[dom]["simplices"]
#             if hlf==1
#                 smplx=get_reflected_index.(smplx)
#             end
#             insert_smplx!(smplx_list, get_rotated_index.(smplx, sctr))
#         end
#
#     elseif sctr>-1 #build sector
#         for smplx in mesh.domains[dom]["simplices"]
#             rsmplx=get_reflected_index.(smplx)
#             insert_smplx!(smplx_list, get_rotated_index.(smplx, sctr))
#             insert_smplx!(smplx_list, get_rotated_index.(rsmplx, sctr))
#         end
#
#     else #build full annulus or just call domain
#         for smplx in mesh.domains[dom]["simplices"]
#             rsmplx=get_reflected_index.(smplx)
#             for sector=0:DOS-1
#                 insert_smplx!(smplx_list, get_rotated_index.(smplx, sector))
#                 insert_smplx!(smplx_list, get_rotated_index.(rsmplx, sector))
#             end
#         end
#     end
#
#
#     domain["simplices"] = smplx_idx
#     return domain
# end
#
#
# import Base.getindex(mesh::Mesh, i::Int)
#
# end


include("sorter.jl")
include("annular_meshes.jl")
include("vtk_write.jl")
##

"""
    points, lines, triangles, tetrahedra, domains=read_msh4(file_name)

Read gmsh's .msh-format version 4.

http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
"""
function read_msh4(file_name)
    lines=Array{UInt32,1}[]
    triangles=Array{UInt32,1}[]
    tetrahedra=Array{UInt32,1}[]
    tag2dom=Dict()
    ent2dom=Array{Dict}(undef,4)
    local points
    domains=Dict()
    open(file_name) do fid
        while !eof(fid)
            line=readline(fid)
            fldname=line[2:end]
            if fldname=="MeshFormat"
                line=readline(fid)
                #if split(line)[1] != "4.1"
                #    break
                #end

            elseif fldname=="PhysicalNames"
                line=readline(fid)
                numPhysicalNames=parse(UInt32,line)
                #Array{String}(undef,numPhysicalNames)
                for idx =1:numPhysicalNames
                    line=readline(fid)
                    dim, tag, dom=split(line)
                    dim=parse(UInt32,dim)
                    dom=dom[2:end-1]
                    tag2dom[tag]=dom
                    domains[dom]=Dict()
                    domains[dom]["dimension"]=dim
                    domains[dom]["simplices"]=[]
                end

            elseif fldname=="Entities"
                line=readline(fid)
                numPoints, numCurves, numSurfaces, numVolumes=parse.(UInt32,split(line))
                for idx=1:4 #inititlaize entity dictionairy
                    ent2dom[idx]=Dict()
                end
                for idx=1:numPoints
                    line=readline(fid)
                    splitted=split(line)
                    entTag=splitted[1]
                    numPhysicalTags=parse(UInt32,splitted[5])
                    physicalTags=splitted[6:5+numPhysicalTags]
                    ent2dom[1][entTag]=[tag2dom[tag] for tag in physicalTags]
                end
                for idx=1:numCurves
                    line=readline(fid)
                    splitted=split(line)
                    entTag=splitted[1]
                    numPhysicalTags=parse(UInt32,splitted[8])
                    physicalTags=splitted[9:8+numPhysicalTags]
                    ent2dom[2][entTag]=[tag2dom[tag] for tag in physicalTags]
                end
                for idx=1:numSurfaces
                    line=readline(fid)
                    splitted=split(line)
                    entTag=splitted[1]
                    numPhysicalTags=parse(UInt32,splitted[8])
                    physicalTags=splitted[9:8+numPhysicalTags]
                    ent2dom[3][entTag]=[tag2dom[tag] for tag in physicalTags]
                end
                for idx=1:numVolumes
                    line=readline(fid)
                    splitted=split(line)
                    entTag=splitted[1]
                    numPhysicalTags=parse(UInt32,splitted[8])
                    physicalTags=splitted[9:8+numPhysicalTags]
                    ent2dom[4][entTag]=[tag2dom[tag] for tag in physicalTags]
                end

            elseif fldname=="Nodes"
                line=readline(fid)
                numEntityBlocks, numNodes,minNodeTag, maxNodeTag =parse.(UInt32,split(line))
                points=Array{Float64}(undef,3,numNodes)
                for idx =1:numEntityBlocks
                    line=readline(fid)
                    entityDim, entityTag, parametric, numNodesInBlock=parse.(UInt32,split(line))
                    #TODO: support parametric
                    nodeTags=Array{UInt32}(undef,numNodesInBlock)
                    for jdx =1:numNodesInBlock
                        line=readline(fid)
                        nodeTags[jdx]=parse(UInt32,line)
                    end
                    for jdx=1:numNodesInBlock
                        line=readline(fid)
                        xyz=parse.(Float64,split(line))
                        points[:,nodeTags[jdx]]=xyz[1:3]
                        #TODO: xyz[4:6] if parametric
                    end
                end
            elseif fldname=="Elements"
                line=readline(fid)
                numEntityBlocks, numElements, minElementTag, maxElementTag = parse.(UInt32,split(line))
                for idx=1:numEntityBlocks
                    line=readline(fid)
                    splitted=split(line)
                    entityDim=parse(UInt32,splitted[1])
                    entityTag=splitted[2]
                    elementType, numElementsInBlock = parse.(UInt32,splitted[3:4])
                    for jdx =1:numElementsInBlock
                        line=readline(fid)
                        nodeTags=parse.(UInt32,split(line))[2:end]
                        if elementType==1
                            append!(lines,[nodeTags]) #TODO: use simplexSorter!!!
                            for dom in ent2dom[entityDim+1][entityTag]
                                append!(domains[dom]["simplices"],length(lines))
                            end
                        elseif elementType==2
                            append!(triangles,[nodeTags]) #TODO: Use simplexSorter!!!
                            for dom in ent2dom[entityDim+1][entityTag]
                                append!(domains[dom]["simplices"],length(triangles))
                            end
                        elseif elementType==4
                            append!(tetrahedra,[nodeTags])
                            for dom in ent2dom[entityDim+1][entityTag]
                                append!(domains[dom]["simplices"],length(tetrahedra))
                            end
                        end
                    end
                end
            end
            #find closing tag and read it
            #while line[2:end]!="End"*fldname
            #    line=readline(fid)
            #end
        end#while
    end#do
    return points, lines, triangles, tetrahedra, domains
end#function


"""
    points, lines, triangles, tetrahedra, domains=read_msh2(file_name)

Read gmsh's .msh-format version 2. This is an old format and wherever possible version 4 should be used.

http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
"""
function read_msh2(file_name)
    lines=Array{UInt32,1}[]
    triangles=Array{UInt32,1}[]
    tetrahedra=Array{UInt32,1}[]
    points=[]
    domains=Dict()
    open(file_name) do f
        field=[]
        field_open=false
        i=1
        #global points, triangles, lines, tetrahedra, domains
        local field_name, domain_map
        while !eof(f)
            line=readline(f)
            if (line[1]=='$')
                if field_open==false
                    field_name=line[2:end]
                    field_open=true
                elseif field_open && (line[2:end]=="End"*field_name)
                   field_open=false
                else
                   println("Error:Field is invalid")
                   break
                end
            end
            if field_open
                if field_name=="MeshFormat"
                    MeshFormat=readline(f)
                elseif field_name=="PhysicalNames"
                    #get number of field entries
                    number_of_entries=parse(Int,readline(f))
                    domain_map=Dict()
                    for i=1:number_of_entries
                        line=split(readline(f))
                        line[3]=strip(line[3],['\"'])
                        domains[line[3]]=Dict("dimension"=>parse(UInt32,line[1])
                        ,"points"=>UInt32[])
                        domain_map[parse(UInt32,line[2])]=line[3]
                    end
                elseif field_name=="Nodes"

                    #get number of field entries
                    number_of_entries=parse(Int,readline(f))
                    points=Array{Float64}(undef,3,number_of_entries)
                        for i=1:number_of_entries
                            line=split(readline(f))
                            points[:,i]=[parse(Float64,line[2]) parse(Float64,line[3]) parse(Float64,line[4])]
                        end
                elseif field_name=="Elements"
                    #get number of field entries
                    number_of_entries=parse(Int,readline(f))
                    for i=1:number_of_entries
                        line=split(readline(f))
                        element_type=parse(Int16,line[2])
                        number_of_tags=parse(Int16,line[3])
                        tags=line[4:4+number_of_tags-1]
                        node_number_list=UInt32[]

                        for node in line[4+number_of_tags:end]
                            append!(node_number_list,parse(UInt32,node))
                        end

                        dom=domain_map[parse(UInt32,tags[1])]
                        for node in node_number_list
                            insert_smplx!(domains[dom]["points"],node)
                        end

                        if element_type==4 #4-node tetrahedron
                        insert_smplx!(tetrahedra,node_number_list)
                        elseif element_type==2 #3-node triangle
                            insert_smplx!(triangles,node_number_list)
                        elseif element_type==1 #2-node line
                            insert_smplx!(lines,node_number_list)
                        #elseif element type==15 #1-node point
                        #    insert_smplx!(      ,node_number_list)
                        else
                            println("error, element type is not defined")
                            return
                        end
                    end
                end
            end
            i+=1
        end
    end
    #triangles=Array{UInt32,1}[]#hot fix for triangle bug
    return points, lines, triangles, tetrahedra, domains
end






##

##

##
"""
    link_triangles_to_tetrahedra!(mesh::Mesh)

Find the tetrahedra that are connected to the triangles in mesh.triangles and store this information in mesh.tri2tet.
"""
function link_triangles_to_tetrahedra!(mesh::Mesh)
    for (idx,tri) in enumerate(mesh.triangles)
        for (idt::UInt32,tet) in enumerate(mesh.tetrahedra)
            if all([i in tet for i=tri])
                mesh.tri2tet[idx]=idt
                break
            end
        end
    end
    return
end







##
"""
    compute_volume!(mesh::Mesh,dom::String)

Compute the size of the domain `dom` and store it in the field `"size"` of the domain definition.
The size will be a length, area, or volume depending on the dimension of the domain.
"""
function compute_volume!(mesh::Mesh,dom::String)
    dim=mesh.domains[dom]["dimension"]
    V=0 #TODO: check whether domains are properly defined
    if dim ==3
        for tet in mesh.tetrahedra[mesh.domains[dom]["simplices"]]
            tet=mesh.points[:,tet]
            V+=abs(LinearAlgebra.det(tet[:,1:3]-tet[:,[4,4,4]]))/6
        end
    #mesh.domains[dom]["volume"]=V
    elseif dim == 2
        for tri in mesh.triangles[mesh.domains[dom]["simplices"]]
            tri=mesh.points[:,tri]
            V+=LinearAlgebra.norm(LinearAlgebra.cross(tri[:,1]-tri[:,3],tri[:,2]-tri[:,3]))/2
        end
    #mesh.domains[dom]["area"]=V
    elseif dim == 1
        for ln in in mesh.lines[mesh.domains[dom]["simplices"]]
            ln=mesh.points[:,ln]
            V+=LinearAlgebra.norm(ln[:,2]-ln[:,1])
        end

    end
    mesh.domains[dom]["size"]=V
end
##
# function get_simplices(mesh::Mesh, dom::String)
#     dim=mesh.domains[dom]["dimension"]
#     if dim==3
#         return mesh.tetrahedra[mesh.domains[dom]["simplices"]]
#     elseif dim==2
#         return mesh.triangles[mesh.domains[dom]["simplices"]]
#     end
# end
##

"""
    tet_idx=find_tetrahedron_containing_point(mesh::Mesh,point)

Find the tetrahedron containing the point `point` and return the index
of the tetrahedron as `tet_idx`. If the point lies on the interface of two or more tetrahedra, the returned `tet_idx` will be the lowest index of these tetrahedra,
i.e. the index depends on the ordering of the tetrahedra in the `mesh.tetrahedra`.
If no tetrahedron in the mesh encloses the point `tet_idx==0`.
"""
function find_tetrahedron_containing_point(mesh::Mesh,point)
    tet_idx=0
    for (idx::UInt32,tet) in enumerate(mesh.tetrahedra)
        tet=mesh.points[:,tet]
        J=tet[:,1:3]-tet[:,[4,4,4]]
        #convert to local barycentric coordinates
        ξ=J\(point-tet[:,4])
        ξ=vcat(ξ,1-sum(ξ))
        if all(0 .<= ξ .<= 1)
            tet_idx=idx
            break
        end
        #TODO impleemnt points on surfaces
    end
    return tet_idx
end

##
function keep!(mesh::Mesh,keys_to_be_kept)
    for key in keys(mesh.domains)
        if !(key in keys_to_be_kept)
            delete!(mesh.domains,key)
        end
    end
end

"""
    collect_lines!(mesh::Mesh)

Populate the list of lines in `mesh.lines`.
"""
function collect_lines!(mesh::Mesh)
    for tet in mesh.tetrahedra
        insert_smplx!(mesh.lines,[tet[1],tet[2]])
        insert_smplx!(mesh.lines,[tet[1],tet[3]])
        insert_smplx!(mesh.lines,[tet[1],tet[4]])
        insert_smplx!(mesh.lines,[tet[2],tet[3]])
        insert_smplx!(mesh.lines,[tet[2],tet[4]])
        insert_smplx!(mesh.lines,[tet[3],tet[4]])
    end
end

"""
    unify!(mesh::Mesh ,new_domain::String,domains...)

Create union of `domains` and name it `new_domain`.
`domains` must share the same dimension and `new_domain`
must be a new name otherwise errors will occur.
"""
function unify!(mesh::Mesh ,new_domain::String,domains...)
    #check new domain is non-existent
    if new_domain in keys(mesh.domains)
        println("ERROR: Domain '"*new_domain*"' already exists!\n Use a different name.")
        return
    end
    #check dimension
    dim=mesh.domains[domains[1]]["dimension"]
    for dom in domains
        if dim != mesh.domains[dom]["dimension"]
            println("Error: Domains are of different dimension!")
            return
        end
    end

    mesh.domains[new_domain]=Dict{String,Any}()
    mesh.domains[new_domain]["dimension"]=dim
    mesh.domains[new_domain]["simplices"]=UInt32[]
    for dom in domains
        #mesh.domains[new_domain]["simplices"]=
        union!(mesh.domains[new_domain]["simplices"],mesh.domains[dom]["simplices"])
    end
    sort!(mesh.domains[new_domain]["simplices"])
end

## legacy functions will be removed in the future
function build_domains!(mesh::Mesh)
    for dom in keys(mesh.domains)
        dim=mesh.domains[dom]["dimension"]
        mesh.domains[dom]["simplices"]=UInt32[]
        if dim == 3
            for (idx::UInt32,tet) in enumerate(mesh.tetrahedra)
                count=0
                for node in tet
                    count+=node in mesh.domains[dom]["points"]
                end
                if count==4
                    append!(mesh.domains[dom]["simplices"],idx)
                end
            end
        elseif dim ==2
            for (idx::UInt32,tri) in enumerate(mesh.triangles)
                count=0
                for node in tri
                    count+=node in mesh.domains[dom]["points"]
                end
                if count==3
                    append!(mesh.domains[dom]["simplices"],idx)
                end
            end
        end
    end
end

function find_tetrahedron_containing_triangle!(mesh::Mesh,dom::String)
    tetmap=Array{UInt32,1}(undef,length(mesh.domains[dom]["simplices"]))
    for (idx,tri) in enumerate(mesh.domains[dom]["simplices"])
        tri=mesh.triangles[tri]
        for (idt::UInt32,tet) in enumerate(mesh.tetrahedra)
            if all([i in tet for i=tri])
                tetmap[idx]=idt
                break
            end
        end
    end
    mesh.domains[dom]["tetmap"]=tetmap
    return nothing
end

""" read a mesh from ansys cfx file"""
function read_ansys(filename)
    lines=Array{UInt32,1}[]
    triangles=Array{UInt32,1}[]
    names=Dict()
    local points,domains,tetrahedra
    #points=[]
    open(filename) do f
        domains=Dict()
        while !eof(f)
            line=readline(f)
            splitted=split(line)
            if line==""
                continue
            end
            if splitted[1]=="(10" && splitted[2]=="(0"
                 #point declaration
                 number_of_points=parse(Int,splitted[4],base=16)
                 points=zeros(3,number_of_points)
                 pnt_idx=1
            elseif splitted[1]=="(10" && splitted[2]!="(0"
                #point field header
                pnt_field=true
                zone_id=splitted[2][2:end]
                first_index=parse(Int,splitted[3],base=16)
                last_index=parse(Int,splitted[4],base=16)
                type=splitted[5]
                #read point field
                if length(splitted)==6
                    ND=parse(UInt32,splitted[6][1:end-2],base=16)
                else
                    ND=-1
                end
                domains[zone_id]=Dict("dimension"=>ND,"points"=>UInt32[])
                for idx=first_index:last_index
                    line=readline(f)
                    splitted=split(line)
                    points[:,idx]=parse.(Float64,splitted)
                    append!(domains[zone_id]["points"],idx)
                end

            elseif splitted[1]=="(12" && splitted[2]=="(0"
                #cell declaration
                number_of_tetrahedra=parse(Int,splitted[4],base=16)
                tetrahedra=Array{Array{UInt32,1}}(undef,number_of_tetrahedra)
                for idx in 1:number_of_tetrahedra
                    tetrahedra[idx]=[]
                end

            elseif splitted[1]=="(13" && splitted[2]=="(0"
                #face declaration
                number_of_faces=parse(Int,splitted[4],base=16)
            elseif splitted[1]=="(13" && splitted[2]!="(0"
                #face field header
                zone_id=splitted[2][2:end]
                first_index=parse(Int,splitted[3],base=16)
                last_index=parse(Int,splitted[4],base=16)
                type=splitted[5]
                element_type=parse(Int,splitted[6][1:end-2],base=16)
                if type=="2"
                    ND=3
                else
                    ND=2
                end
                domains[zone_id]=Dict("dimension"=>ND,"points"=>UInt32[])

                #read face field
                for idx =first_index:last_index
                    line=readline(f)
                    splitted=split(line)
                    left=parse(UInt32,splitted[end-1],base=16)
                    right=parse(UInt32,splitted[end],base=16)
                    tri=parse.(UInt32,splitted[1:3],base=16)
                    for pnt in tri
                        #Hier können die Domain simplices abgefragt werden
                        if left!=0
                            insert_smplx!(tetrahedra[left],pnt)
                        end
                        if right!=0
                            insert_smplx!(tetrahedra[right],pnt)
                        end
                        insert_smplx!(domains[zone_id]["points"],pnt)
                    end
                    if left==0 || right==0
                        insert_smplx!(triangles,tri)
                        #Das muss Rückverfolgt werden um die Dreiecke zu identifizieren
                    end
                end
            elseif splitted[1]=="(45"
                names[splitted[2][2:end]]=splitted[end][1:end-4]
            end
        end
    end
    for (key,value) in names
        if key in keys(domains)
            domains[value]=pop!(domains,key)
        else
            #TODO: warning message
        end

    end
    #points, lines, triangles, tetrahedra, domains, names
    return Mesh(filename, points, lines, triangles, tetrahedra, domains, filename,1)
end


end#module
