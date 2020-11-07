# open(filename) do fid
#     while !eof(fid)
#         line = readline(fid)
#     end
# end
# ##
# ##
# fid = open("shape_sense_v1_lin.vtu")
# ##
# line = readline(fid)
# splitted=split(strip(line),('>','<'))
# if length(splitted[2])!=0
#     println(splitted)
# end
# close(fid)
# ##
# using LightXML
#
# xdoc = parse_file("shape_sense_v1_lin.vtu")
#
# xroot = root(xdoc)
#
# ces = get_elements_by_tagname(xroot, "UnstructuredGrid")

##
using EzXML
using XMLDict

function read_data_from_vtu(filename)
    xdoc = readxml(filename)
    D=xml_dict(xdoc)
    data=D["VTKFile"]["UnstructuredGrid"]["Piece"]["PointData"]["DataArray"]#[1][""]
    Data=Dict()
    for idx in 1:length(data)
        Data[data[idx][:Name]]=parse.(Float64,split(data[idx][""]))
    end
    return Data
end
##
using WavesAndEigenvalues.Helmholtz
##
fname="./examples/tutorials/Rijke_hot_lin.vtu"
fname2="./examples/tutorials/Rijke_hot_fd_lin.vtu"
data=read_data_from_vtu(fname)
data2=read_data_from_vtu(fname2)

##
for coo in ("x","y","z")
    mode="real DA $coo[272.27 + 1.53im]Hz"
    #mode="real DA $coo[520.73 + 0.0im]Hz"
    devi=data[mode].-data2[mode]
    data["devi $coo"]=devi
end
mesh=Mesh("./examples/tutorials/Rijke_mm.msh",scale=0.001)


vtk_write("Rijke_hot_devi",mesh, data)

##
xdoc = readxml("./examples/shape/shape_sense_v1_lin.vtu")

D=xml_dict(xdoc)


data=D["VTKFile"]["UnstructuredGrid"]["Piece"]["PointData"]["DataArray"]#[1][""]
D["VTKFile"]["UnstructuredGrid"]["Piece"]["Cells"]["DataArray"]
##
points=D["VTKFile"]["UnstructuredGrid"]["Piece"]["Points"]["DataArray"][""]

points=parse.(Float64,split(points))
points=reshape(points,3,length(points)÷3)
##
cells=D["VTKFile"]["UnstructuredGrid"]["Piece"]["Cells"]["DataArray"]#[""]

cells=parse.(Int,split(cells[1][""]))
cells=reshape(cells,4,length(cells)÷4)

##
data=D["VTKFile"]["UnstructuredGrid"]["Piece"]["PointData"]["DataArray"]#[1][""]
Data=Dict()
for idx in 1:length(data)
    Data[data[idx][:Name]]=parse.(Float64,split(data[idx][""]))
end
