"""
Module containing functionality to read and process tetrahedral meshes in gmsh or nastran format.
"""
module Meshutils
import LinearAlgebra, ProgressMeter
include("Meshutils_exports.jl")

"""
Simple struct element to store additional information on symmetry for rotational symmetric meshes.

#Fields
- `DOS::Int64`: degree of symmetry
- `naxis::Int64`: number of grid points lying on the center axis
- `nxbloch::Int64`: number of grid points lying on the Bloch plane but not on the center axis
- `nbody::Int64`: number of interior points in half cell
- `shiftbody::Int64`: difference from body point to reflected body point
- `nxsymmetry::Int64`: number of grid points on symmetry plane (without center axis)
- `nxsector::Int64`: number of grid points belonging to a unit cell (without cenetraxis, Bloch and Bloch image plane)
- `naxis_ln::Int64`: number of line segments lying on the center axis
- `nxbloch_ln::Int64`: number of line segments belonging to a unit cell (without cenetraxis, Bloch and Bloch image plane)
- `nxsector_ln::Int64`: number of line segments belonging to a unit cell (without cenetraxis, Bloch and Bloch image plane)
- `nxsector_tri::Int64`: number of surface triangles belonging to a unit cell (without Bloch and Bloch image plane)
- `nxsector_tet::Int64`: number of tetrahedra of a unit cell
- `n`: unit axis vector of the center axis
- `pnt: foot point of the center axis
- `unit::Bool`: true if mesh represents only a unit cell
"""
struct SymInfo
    DOS::Int64 #degree of symmetry
    naxis::Int64 #number of gridpoints lying on the centeraxis
    nxbloch::Int64
    nbody::Int64
    shiftbody::Int64
    nxsymmetry::Int64
    nxsector::Int64
    naxis_ln::Int64
    nxbloch_ln::Int64
    nxsector_ln::Int64
    nxsector_tri::Int64
    nxsector_tet::Int64
    n
    pnt
    unit::Bool
end


"""
Definition of the mesh type

# Fields
- `name::String`: the name of the mesh.
- `points::Array`: 3×N array containing the coordinates of the N points defining the mesh.
- `lines::List`: List of simplices defining the edges of the mesh
- `triangles::List`: List of simplices defining the surface triangles of the mesh
- `int_triangles::List`: List of simplices defining the interior triangles of the mesh
- `tetrahedra::List`: List of simplices defining the tetrahedra of the mesh
- `domains::Dict`: Dictionary defining the domains of the mesh. See comments below.
- `file::String`: path to the file containing the mesh.
- `tri2tet::Array`: Array of length `length(tetrahedra)` containing the indices of the connected tetrahedra.
- `dos`: special field meant to contain symmetry information of highly symmetric meshes.

# Notes
The meshes are supposed to be tetrahedral. All simplices (lines, triangles, and tetrahedra) are stored as lists of simplices.
Simplices are lists of integers containing the indices of the points (i.e. the column number in the `points` array) forming the simplex.
This means a line is a two-entry list, a triangle a three-entry list, and a tetrahedron a four-entry list.
For convenience certain entities of the mesh can be further defined in the `domains` dictionary. Each key defines a domain and maps to another dictionary.
This second-level dictionary contains at least two keys: `"dimension"` mapping to the dimension of the specified domain (1,2, or 3) and
`"simplices"` containing a list of integers mapping into the respective simplex lists. More keys may be added to the dictionary to
define additional and/or custom information on the domain. For instance the `compute_size!` function adds an entry with the domain size.
"""
struct Mesh
    name
    points
    lines
    triangles
    int_triangles
    tetrahedra
    domains
    file
    tri2tet
    dos
end
## constructor
"""
    mesh=Mesh(file_name::String; scale=1, inttris=false)

read a tetrahedral mesh from gmsh or nastran file into `mesh`.
The optional scaling factor `scale` may be used to scale the units of the mesh.
The optional parameter `inttris` toggles whether triangles are sortet into two
seperate lists for surface and internal triangles.
"""
function Mesh(file_name::String;scale=1,inttris=false)
    #mesh constructor
    ext=split(file_name,".")
    if ext[end]=="msh"
        #gmsh mesh file
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
    elseif ext[end]=="nas" || ext[end]=="bdf"
        points, lines, triangles, tetrahedra, domains =read_nastran(file_name)
        # if isfile(ext[1]*"_surfaces.txt")
        #     relabel_nas_domains!(domains,ext[1]*"_surfaces.txt")
        # end
    else
        println("mesh type not supported")
        return nothing
    end

    #sometimes elements are multiply defined, make them unique
    #TODO: This can be done more memory efficient by reading the msh-file twice
    ulines=Array{UInt32,1}[]
    utriangles=Array{UInt32,1}[]
    utetrahedra=Array{UInt32,1}[]
    uinttriangles=Array{UInt32,1}[]

    for ln in lines
        insert_smplx!(ulines,ln)
    end

    for tet in tetrahedra
        insert_smplx!(utetrahedra,tet)
    end

    if inttris
        utriangles,uinttriangles=assemble_triangles(utetrahedra)
    else
        for tri in triangles
            insert_smplx!(utriangles,tri)
        end
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

    if inttris
        tri2tet=[zeros(UInt32,length(utriangles)),zeros(UInt32,length(uinttriangles),2)]
    else
        tri2tet=zeros(UInt32,length(utriangles))
        if length(tri2tet)!=0
            tri2tet[1]=0xffffffff #magic sentinel value to check whether list is linked
        end
    end
    Mesh(file_name, points*scale, ulines, utriangles, uinttriangles, utetrahedra, domains, file_name,tri2tet,1)
end
#deprecated constructor
function Mesh(file, points, lines, triangles, tetrahedra, domains, file_name, tri2tet, dos)
    return Mesh(file, points, lines, triangles, [], tetrahedra, domains, file_name, tri2tet, dos)
end

#show mesh
import Base.string
function string(mesh::Mesh)
    if isempty(mesh.int_triangles)
        txt="""mesh: $(mesh.file)
        #################
        points:     $(size(mesh.points)[2])
        lines:      $(length(mesh.lines))
        triangles:  $(length(mesh.triangles))
        tetrahedra: $(length(mesh.tetrahedra))
        #################
        domains: """
    else
        txt="""mesh: $(mesh.file)
        #################
        points:            $(size(mesh.points)[2])
        lines:             $(length(mesh.lines))
        triangles (surf):  $(length(mesh.triangles))
        triangles (int):   $(length(mesh.int_triangles))
        tetrahedra:        $(length(mesh.tetrahedra))
        #################
        domains: """
    end
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


include("./Mesh/sorter.jl")
include("./Mesh/annular_meshes.jl")
include("./Mesh/vtk_write.jl")
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
                    dom=String(dom[2:end-1])
                    tag2dom[tag]=dom
                    domains[dom]=Dict()
                    domains[dom]["dimension"]=dim
                    domains[dom]["simplices"]=[]
                end

            elseif fldname=="Entities"
                line=readline(fid)
                numPoints, numCurves, numSurfaces, numVolumes=parse.(UInt32,split(line))
                for idx=1:4 #inititlaize entity dictionary
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

include("./Mesh/read_nastran.jl")




##

##

##
"""
    link_triangles_to_tetrahedra!(mesh::Mesh)

Find the tetrahedra that are connected to the triangles in mesh.triangles (and mesh.inttriangles) and store this information in mesh.tri2tet.
"""
function link_triangles_to_tetrahedra!(mesh::Mesh)
    if !isempty(mesh.int_triangles)
        #IDEA: move this to constructor
        for (idt,tet) in enumerate(mesh.tetrahedra)
            for tri in (tet[[1,2,3]],tet[[1,2,4]],tet[[1,3,4]],tet[[2,3,4]])
                idx,ins=sort_smplx(mesh.triangles,tri)
                if !ins
                    mesh.tri2tet[1][idx]=idt
                else
                    idx,ins=sort_smplx(mesh.int_triangles,tri)
                    if !ins
                        if mesh.tri2tet[2][idx,1]==0
                            mesh.tri2tet[2][idx,1]=idt
                        else
                            mesh.tri2tet[2][idx,2]=idt
                        end
                    end
                end
            end
        end

    else
        for (idt,tet) in enumerate(mesh.tetrahedra)
            for tri in (tet[[1,2,3]],tet[[1,2,4]],tet[[1,3,4]],tet[[2,3,4]])
                idx,ins=sort_smplx(mesh.triangles,tri)
                if !ins
                    mesh.tri2tet[idx]=idt
                end
            end
        end
    end
    return nothing
end
##
function assemble_triangles(tetrahedra)
    int_triangles=Array{UInt32,1}[]
    ext_triangles=Array{UInt32,1}[]
    for tet in tetrahedra#loop over all tetrahedra
        for tri_idcs in [[1,2,3],[1,2,4],[1,3,4],[2,3,4]] #loop over all 4 triangles of the tetrahedron
            tri=tet[tri_idcs]
            idx,notinside=sort_smplx(ext_triangles, tri)
            if notinside #triangle occurs for the first time
                #IDEA: there is no sanity check. If tetrahedra are ill defined a
                #triangle might show up three times. That would be fatal.
                #may implement some low-level-sanity checks ensuring data
                #integrity

                insert!(ext_triangles,idx,tri)
            else #triangle occurs for the second time, it must be an interior triangle.
                deleteat!(ext_triangles,idx)
                insert_smplx!(int_triangles,tri)
            end
        end
    end
    return ext_triangles,int_triangles
end
##
"""
    new_mesh=octosplit(mesh::Mesh)

Subdivide each tetrahedron in the mesh `mesh` into 8 tetrahedra by splitting
each edge at its center.

# Notes
The algorithm introduces `length(mesh.lines)` new vertices. This yields a finer
mesh featuring `8*size(mesh,2)` tetrahedra, `4*length(mesh.triangles)`
triangles, and `2*length(mesh.lines)` lines. From the  3 possible subdivision of
a tetrahedron the algorithm automatically chooses the one that minimizes the
edge lengths of the new tetrahedra. The point labeling of `new_mesh` is
consistent with the point labeling in `mesh`, i.e., the first
`size(mesh.points,2)` points in `new_mesh.points` are identical to the points
in `mesh.points`. Hence, `mesh` and `new_mesh` form a hierachy of meshes.
"""
function octosplit(mesh::Mesh)
    collect_lines!(mesh)
    N_points=size(mesh.points,2)
    N_lines=length(mesh.lines)
    N_tet=length(mesh.tetrahedra)
    N_tri=length(mesh.triangles)
    ##
    points=Array{Float64,2}(undef,3,N_points+N_lines)
    ## extend point list
    points[:,1:N_points]=mesh.points
    for (idx,ln) in enumerate(mesh.lines)
        points[:,N_points+idx]= sum(mesh.points[:,ln],dims=2)./2
    end
    #create new tetrahedra by splitting the old ones
    tetrahedra=Array{UInt32}[]
    for tet in mesh.tetrahedra
        A,B,C,D=tet #old vertices
        #indices of the center points of the lines also become vertices
        AB=N_points+find_smplx(mesh.lines,[A,B])
        AC=N_points+find_smplx(mesh.lines,[A,C])
        AD=N_points+find_smplx(mesh.lines,[A,D])
        BC=N_points+find_smplx(mesh.lines,[B,C])
        BD=N_points+find_smplx(mesh.lines,[B,D])
        CD=N_points+find_smplx(mesh.lines,[C,D])

        #outer tetrahedra always present...
        insert_smplx!(tetrahedra,[A, AB, AC, AD])
        insert_smplx!(tetrahedra,[B, AB, BC, BD])
        insert_smplx!(tetrahedra,[C, AC, BC, CD])
        insert_smplx!(tetrahedra,[D, AD, BD, CD])

        # split remaining octet based on minimum distance
        AB_CD=LinearAlgebra.norm(points[:,AB].-points[:,CD])
        AC_BD=LinearAlgebra.norm(points[:,AC].-points[:,BD])
        AD_BC=LinearAlgebra.norm(points[:,AD].-points[:,BC])

        if AB_CD<=AC_BD && AB_CD<=AD_BC
            insert_smplx!(tetrahedra,[AB,CD,AC,AD])
            insert_smplx!(tetrahedra,[AB,CD,AD,BD])
            insert_smplx!(tetrahedra,[AB,CD,BD,BC])
            insert_smplx!(tetrahedra,[AB,CD,BC,AC])
        elseif AC_BD<=AB_CD && AC_BD<=AD_BC
            insert_smplx!(tetrahedra,[AC,BD,AB,AD])
            insert_smplx!(tetrahedra,[AC,BD,AD,CD])
            insert_smplx!(tetrahedra,[AC,BD,CD,BC])
            insert_smplx!(tetrahedra,[AC,BD,BC,AB])
        elseif AD_BC<=AC_BD && AD_BC<=AB_CD
            insert_smplx!(tetrahedra,[AD,BC,AC,CD])
            insert_smplx!(tetrahedra,[AD,BC,CD,BD])
            insert_smplx!(tetrahedra,[AD,BC,BD,AB])
            insert_smplx!(tetrahedra,[AD,BC,AB,AC])
        end
    end

    #split the old surface triangles
    triangles=Array{UInt32}[]
    for tri in mesh.triangles
        A,B,C=tri
        AB=N_points+find_smplx(mesh.lines,[A,B])
        AC=N_points+find_smplx(mesh.lines,[A,C])
        BC=N_points+find_smplx(mesh.lines,[B,C])
        insert_smplx!(triangles,[A,AB,AC])
        insert_smplx!(triangles,[B,AB,BC])
        insert_smplx!(triangles,[C,AC,BC])
        insert_smplx!(triangles,[AB,AC,BC])
    end

    #split the old interior triangles
    #TODO:Testing
    inttriangles=Array{UInt32}[]
    # for tri in mesh.int_triangles
    #     A,B,C=tri
    #     AB=N_points+find_smplx(mesh.lines,[A,B])
    #     AC=N_points+find_smplx(mesh.lines,[A,C])
    #     BC=N_points+find_smplx(mesh.lines,[B,C])
    #     insert_smplx!(inttriangles,[A,AB,AC])
    #     insert_smplx!(inttriangles,[B,AB,BC])
    #     insert_smplx!(inttriangles,[C,AC,BC])
    #     insert_smplx!(inttriangles,[AB,AC,BC])
    # end

    ## relabel tetrahedra
    tet_labels=Array{UInt32}(undef,N_tet,8)
    for (idx,tet) in enumerate(mesh.tetrahedra)
        A,B,C,D=tet #old vertices
        #indices of the center points of the lines also become vertices
        AB=N_points+find_smplx(mesh.lines,[A,B])
        AC=N_points+find_smplx(mesh.lines,[A,C])
        AD=N_points+find_smplx(mesh.lines,[A,D])
        BC=N_points+find_smplx(mesh.lines,[B,C])
        BD=N_points+find_smplx(mesh.lines,[B,D])
        CD=N_points+find_smplx(mesh.lines,[C,D])

        #outer tetrahedra always present...
        tet_labels[idx,1]=find_smplx(tetrahedra,[A, AB, AC, AD])
        tet_labels[idx,2]=find_smplx(tetrahedra,[B, AB, BC, BD])
        tet_labels[idx,3]=find_smplx(tetrahedra,[C, AC, BC, CD])
        tet_labels[idx,4]=find_smplx(tetrahedra,[D, AD, BD, CD])

        # split remaining octet based on minimum distance
        AB_CD=LinearAlgebra.norm(points[:,AB].-points[:,CD])
        AC_BD=LinearAlgebra.norm(points[:,AC].-points[:,BD])
        AD_BC=LinearAlgebra.norm(points[:,AD].-points[:,BC])

        if AB_CD<=AC_BD && AB_CD<=AD_BC
            tet_labels[idx,5]=find_smplx(tetrahedra,[AB,CD,AC,AD])
            tet_labels[idx,6]=find_smplx(tetrahedra,[AB,CD,AD,BD])
            tet_labels[idx,7]=find_smplx(tetrahedra,[AB,CD,BD,BC])
            tet_labels[idx,8]=find_smplx(tetrahedra,[AB,CD,BC,AC])
        elseif AC_BD<=AB_CD && AC_BD<=AD_BC
            tet_labels[idx,5]=find_smplx(tetrahedra,[AC,BD,AB,AD])
            tet_labels[idx,6]=find_smplx(tetrahedra,[AC,BD,AD,CD])
            tet_labels[idx,7]=find_smplx(tetrahedra,[AC,BD,CD,BC])
            tet_labels[idx,8]=find_smplx(tetrahedra,[AC,BD,BC,AB])
        elseif AD_BC<=AC_BD && AD_BC<=AB_CD
            tet_labels[idx,5]=find_smplx(tetrahedra,[AD,BC,AC,CD])
            tet_labels[idx,6]=find_smplx(tetrahedra,[AD,BC,CD,BD])
            tet_labels[idx,7]=find_smplx(tetrahedra,[AD,BC,BD,AB])
            tet_labels[idx,8]=find_smplx(tetrahedra,[AD,BC,AB,AC])
        end
    end

    #relabel triangles
    tri_labels=Array{UInt32}(undef,N_tri,4)
    for (idx,tri) in enumerate(mesh.triangles)
        A,B,C=tri
        AB=N_points+find_smplx(mesh.lines,[A,B])
        AC=N_points+find_smplx(mesh.lines,[A,C])
        BC=N_points+find_smplx(mesh.lines,[B,C])
        tri_labels[idx,1]=find_smplx(triangles,[A,AB,AC])
        tri_labels[idx,2]=find_smplx(triangles,[B,AB,BC])
        tri_labels[idx,3]=find_smplx(triangles,[C,AC,BC])
        tri_labels[idx,4]=find_smplx(triangles,[AB,AC,BC])
    end

    #reassemble domains
    domains=Dict()
    for (dom,domain) in mesh.domains
        domains[dom]=Dict()
        dim=domain["dimension"]
        domains[dom]["dimension"]=dim
        #TODO: detect whether there optional field like "size"
        simplices=[]
        for smplx_idx in domain["simplices"]
            if dim==3
                append!(simplices,tet_labels[smplx_idx,:])
            elseif dim==2
                append!(simplices,tri_labels[smplx_idx,:])
            end
        end
        domains[dom]["simplices"]=sort(simplices)
    end

    tri2tet=zeros(UInt32,length(triangles))
    if length(tri2tet)>0
        tri2tet[1]=0xffffffff
    end #magic sentinel value to check whether list is linked
    return Mesh(mesh.file, points, [], triangles, inttriangles, tetrahedra, domains, mesh.file,tri2tet,1)
end


##
"""
    compute_size!(mesh::Mesh,dom::String)

Compute the size of the domain `dom` and store it in the field `"size"` of the domain definition.
The size will be a length, area, or volume depending on the dimension of the domain.
"""
function compute_size!(mesh::Mesh,dom::String)
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


"""
    surface_points, tri_mask, tet_mask = get_surface_points(mesh::Mesh,output=true)

Get a list `surface_points` of all point indices that are on the surface of the
mesh `mesh`. The corresponding lists `tri_mask and `tet_mask` contain lists of
indices of all triangles and tetrahedra that are connected to the surface
point. The optional parameter `output` toggles whether a progressbar should be
shown or not.
"""
function get_surface_points(mesh::Mesh,output::Bool=true)
    #TODO put return values as fields into Mesh
    #step #1 find all points on the surface
    surface_points=[]
    #intialize progressbar
    if output
        p = ProgressMeter.Progress(length(mesh.triangles),desc="Find points #1... ", dt=1,
             barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=20)
    end
    for tri in mesh.triangles
        for pnt in tri
            insert_smplx!(surface_points,pnt)
        end
        #update progressbar
        if output
            ProgressMeter.next!(p)
        end
    end

    #step 2 link points to simplices
    #tri_mask=Array{Int64,1}[]
    #tet_mask=Array{Int64,1}[]
    tri_mask=[[] for idx = 1:length(surface_points)]
    tet_mask=[[] for idx = 1:length(surface_points)]
    #intialize progressbar
    if output
        p = ProgressMeter.Progress(length(mesh.triangles),desc="Link to triangles... ", dt=1,
             barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=20)
    end

    for (tri_idx,tri) in enumerate(mesh.triangles)
        for pnt_idx in tri
            idx=find_smplx(surface_points,pnt_idx)
            if idx>0
                append!(tri_mask[idx],tri_idx)
            end
        end
        #update progressbar
        if output
            ProgressMeter.next!(p)
        end
    end
    if output
        p = ProgressMeter.Progress(length(mesh.tetrahedra),desc="Link to tetrahedra... ", dt=1,
             barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=20)
    end
    for (tet_idx,tet) in enumerate(mesh.tetrahedra)
        for pnt_idx in tet
            idx=find_smplx(surface_points,pnt_idx)
            if idx>0
                append!(tet_mask[idx],tet_idx)
            end
        end
        #update progressbar
        if output
            ProgressMeter.next!(p)
        end
    end

    if mesh.dos!=1 && mesh.dos.unit
        #blochify surface points
        for idx in surface_points
            bidx=idx -mesh.dos.naxis
            #reference surface
            if 0 < bidx <= mesh.dos.nxbloch
                append!(tri_mask[idx],tri_mask[end-mesh.dos.nxbloch+bidx])
                append!(tet_mask[idx],tet_mask[end-mesh.dos.nxbloch+bidx])
                #unique!(tri_mask)
                #unique!(tet_mask)
                unique!(tri_mask[idx])
                unique!(tet_mask[idx])
                #symmetrize
                append!(tri_mask[end-mesh.dos.nxbloch+bidx],tri_mask[idx])
                append!(tet_mask[end-mesh.dos.nxbloch+bidx],tet_mask[idx])
                unique!(tri_mask[end-mesh.dos.nxbloch+bidx])
                unique!(tet_mask[end-mesh.dos.nxbloch+bidx])
            end
        end
    end

    return surface_points, tri_mask, tet_mask
end
# function get_surface_points(mesh::Mesh,output::Bool=true)
#     #TODO put return values as fields into Mesh
#     #step #1 find all points on the surface
#     surface_points=[]
#     #intialize progressbar
#     if output
#         p = ProgressMeter.Progress(length(mesh.triangles),desc="Find points #1... ", dt=1,
#              barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
#              barlen=20)
#     end
#     for tri in mesh.triangles
#         for pnt in tri
#             insert_smplx!(surface_points,pnt)
#         end
#         #update progressbar
#         if output
#             ProgressMeter.next!(p)
#         end
#     end
#
#     #step 2 link points to simplices
#     tri_mask=Array{Int64,1}[]
#     tet_mask=Array{Int64,1}[]
#     #intialize progressbar
#     if output
#         p = ProgressMeter.Progress(length(surface_points),desc="Link to simplices... ", dt=1,
#              barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
#              barlen=20)
#     end
#     for (idx,pnt_idx) in enumerate(surface_points)
#         mask=Int64[]
#         for (tri_idx,tri) in enumerate(mesh.triangles)
#             if pnt_idx in tri
#                 append!(mask,tri_idx)
#             end
#         end
#         append!(tri_mask,[mask])
#         mask=Int64[]
#         for (tet_idx,tet) in enumerate(mesh.tetrahedra)
#             if pnt_idx in tet
#                 append!(mask,tet_idx)
#             end
#         end
#         append!(tet_mask, [mask])
#         #update progressbar
#         if output
#             ProgressMeter.next!(p)
#         end
#     end
#     return surface_points, tri_mask, tet_mask
# end

"""
    normal_vectors=get_normal_vectors(mesh::Mesh,output::Bool=true)

Compute a 3×`length(mesh.triangles)` Array containing the outward pointing
normal vectors of each of the surface triangles of the mesh `mesh`. The vectors
are not normalised but their norm is twice the area of the corresponding
triangle. The optional parameter `output` toggles whether a progressbar should be
shown or not.
"""
function get_normal_vectors(mesh::Mesh,output::Bool=true)
    #TODO put return values as fields into Mesh
    #initialize progressbar
    if output
        p = ProgressMeter.Progress(length(mesh.triangles),desc="Find points #1... ", dt=1,
             barglyphs=ProgressMeter.BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
             barlen=20)
    end
    normal_vectors=Array{Float64}(undef,3,length(mesh.triangles))
    #Identify points A,B,C, and D that form the tetrahedron connected to the
    #surface triangle. D is not part of the surface but is inside the body
    #therefore its positio can be used to find determine the outward direction.
    if length(mesh.tri2tet)!=0 && mesh.tri2tet[1]==0xffffffff #TODO: Initialize as empty instead with sentinel value
        link_triangles_to_tetrahedra!(mesh)
    end
    for (idx,tri) in enumerate(mesh.triangles)
        A,B,C=tri
        tet=mesh.tri2tet[idx]
        tet=mesh.tetrahedra[tet]
        D=0
        for idx in tet
            if !(idx in tri)
                D=idx
                break
            end
        end
        A=mesh.points[:,A]
        B=mesh.points[:,B]
        C=mesh.points[:,C]
        D=mesh.points[:,D]
        N=LinearAlgebra.cross(A.-C,B.-C) #compute normal vector (length is twice the area of teh surface triangle)
        N*=sign(LinearAlgebra.dot(N,C.-D))# make normal outward
        normal_vectors[:,idx]=N
    end
    #update progressbar
    if output
        ProgressMeter.next!(p)
    end
    return normal_vectors
end

"""
    generate_field(mesh::Mesh,func;el_type=0)

Generate field from function `func`for mesh `mesh`. The element type is either
`el_type=0` for field values associated with the mesh tetrahedra or `el_type=1`
for field values associated with the mesh vertices. The function must accept
three input arguments corresponding to the three space dimensions.
"""
function generate_field(mesh::Mesh,func;order::Symbol=:const)
    #TODO: implement el_type=2
    if order==:const
        field=zeros(Float64,length(mesh.tetrahedra))#*347.0
        for idx=1:length(mesh.tetrahedra)
            cntr=sum(mesh.points[:,mesh.tetrahedra[idx]],dims=2)/4 #centerpoint of a tetrahedron is arithmetic mean of the vertices
            field[idx]=func(cntr...)
        end
    elseif order==:lin
        field=zeros(Float64,size(mesh.points,2))
        for idx in 1:size(mesh.points,2)
            pnt=mesh.points[:,idx]
            field[idx]=func(pnt...)
        end
    else
        println("Error: Can't generate field order $order not supported!")
        return nothing
    end
    return field
end
##legacycode
# function generate_field(mesh::Mesh,func;el_type=0)
#     #TODO: implement el_type=2
#     if el_type==0
#         field=zeros(Float64,length(mesh.tetrahedra))#*347.0
#         for idx=1:length(mesh.tetrahedra)
#             cntr=sum(mesh.points[:,mesh.tetrahedra[idx]],dims=2)/4 #centerpoint of a tetrahedron is arithmetic mean of the vertices
#             field[idx]=func(cntr...)
#         end
#     else
#         field=zeros(Float64,size(mesh.points,2))
#         for idx in 1:size(mesh.points,2)
#             pnt=mesh.points[:,idx]
#             field[idx]=func(pnt...)
#         end
#     end
#     return field
# end
"""
    data,surf_keys,vol_keys=color_domains(mesh::Mesh,domains=[])

Create data fields containing an integer number corresponding to the local
domain.

# Arguments
- `mesh::Mesh`: mesh to be colored
- `domains::List`: (optional) list of domains that are to be colored. If empty  all domains will be colored.

# Returns
- `data::Dict`: Dictionary containing a key for each colored domain plus magic keys `"__all_surfaces__"` and `"__all_volumes__"` holding all colors in one field. These magic fields may not work correctly if the domains are not disjoint.
- `surf_keys::Dict`: Dictionary mapping the surface domain names to their indeces.
- `vol_keys::Dict`: Dictionary mapping the volume domain names to their indeces.

# Notes
The `data`variable is designed to flawlessly work with `vtk_write` for
visualization with paraview. Check the "Interpret Values as Categories" box in
paraview's color map editor for optimal visualization.

See also: [`vtk_write`](@ref)
"""
function color_domains(mesh::Mesh,domains=[])
    number_of_triangles=length(mesh.triangles)
    number_of_tetrahedra=length(mesh.tetrahedra)
    #IDEA: check number of domains and choose datatype accordingly (UInt8 might be possible)
    tri_color=zeros(UInt16,number_of_triangles)
    tet_color=zeros(UInt16,number_of_tetrahedra)
    data=Dict()
    surf_keys=Dict()
    vol_keys=Dict()
    if length(domains)==0
        domains=sort([key for key in keys(mesh.domains)])
    end
    surf_idx=0
    vol_idx=0
    for key in domains
        if key in keys(mesh.domains)
            dom=mesh.domains[key]
        else
            println("Warning: No domain named '$key' in mesh.")
            continue
        end
        smplcs=dom["simplices"]
        dim=dom["dimension"]
        if dim==2
            surf_idx+=1

            if any(tri_color[smplcs].!=0)
                println("domain $key is overlapping")
            end
            tri_color[smplcs].=surf_idx
            data[key]=zeros(Int64,number_of_triangles)
            data[key][smplcs].=surf_idx
            surf_keys[key]=surf_idx
            println("surface:$key-->$surf_idx")
        elseif dim==3
            vol_idx+=1
            if any(tet_color[smplcs].!=0)
                println("domain $key is overlapping")
            end
            tet_color[smplcs].=vol_idx
            data[key]=zeros(Int64,number_of_tetrahedra)
            data[key][smplcs].=vol_idx
            vol_keys[key]=vol_idx
            println("volume:$key-->$vol_idx")
        end
    end

    data["__all_surfaces__"]=tri_color
    data["__all_volumes__"]=tet_color

 return data,surf_keys,vol_keys
end

#################################################
#! legacy functions will be removed in the future
#!################################################
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
