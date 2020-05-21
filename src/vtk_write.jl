#functionality in this file implements writing .vtu files
#complient tohttps://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
##
function vtk_write1(file_name::String,mesh::Mesh,data=[])
   open(file_name*".vtu","w") do fid
      write(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" >\n")#byte_order=\"BigEndian\"
      write(fid,"\t<UnstructuredGrid>\n")
      write(fid,"\t\t<Piece NumberOfPoints=\"$(size(mesh.points,2))\" NumberOfCells=\"$(length(mesh.tetrahedra))\">\n")
      write(fid,"\t\t\t<Points>\n")
      write(fid,"\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")
      for pnt in mesh.points
         write(fid,string(pnt)*" ")
      end
      write(fid,"</DataArray>")
      write(fid,"\t\t\t</Points>\n")
      write(fid,"\t\t\t<Cells>\n")
      write(fid,"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">")
      for tet in mesh.tetrahedra
         for pnt in tet
            write(fid,string(Int32(pnt-1))*" ")
         end
      end
      write(fid,"</DataArray>\n")
      write(fid,"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">")
      off=0
      for tet in mesh.tetrahedra
            off+=length(tet)
            write(fid,string(off)*" ")
      end
      write(fid,"</DataArray>\n")
      write(fid,"\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">")
      for tet in mesh.tetrahedra
         for pnt in tet
            write(fid,"10 ")
         end
      end
      write(fid,"</DataArray>\n")
      write(fid,"\t\t\t</Cells>\n")

      write(fid,"\t\t\t<PointData>\n")
      for dataname in keys(data)
         write(fid,"\t\t\t\t<DataArray type=\"Float64\" Name=\"$dataname\" format=\"ascii\">")
         for datum in data[dataname]
               write(fid,string(datum)*" ")

         end
         write(fid,"\n\t\t\t\t</DataArray>\n")
      end
      write(fid,"\t\t\t</PointData>\n")
      write(fid,"\t\t</Piece>\n")
      write(fid,"\t</UnstructuredGrid>\n")
      write(fid,"</VTKFile>\n")
   end
   return nothing
end
##
function vtk_write0(file_name::String,mesh::Mesh,data=[])
   open(file_name*".vtu","w") do fid
      write(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" >\n")#byte_order=\"BigEndian\"
      write(fid,"\t<UnstructuredGrid>\n")
      write(fid,"\t\t<Piece NumberOfPoints=\"$(size(mesh.points,2))\" NumberOfCells=\"$(length(mesh.tetrahedra))\">\n")
      write(fid,"\t\t\t<Points>\n")
      write(fid,"\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")
      for pnt in mesh.points
         write(fid,string(pnt)*" ")
      end
      write(fid,"</DataArray>")
      write(fid,"\t\t\t</Points>\n")
      write(fid,"\t\t\t<Cells>\n")
      write(fid,"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">")
      for tet in mesh.tetrahedra
         for pnt in tet
            write(fid,string(Int32(pnt-1))*" ")
         end
      end
      write(fid,"</DataArray>\n")
      write(fid,"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">")
      off=0
      for tet in mesh.tetrahedra
            off+=length(tet)
            write(fid,string(off)*" ")
      end
      write(fid,"</DataArray>\n")
      write(fid,"\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">")
      for tet in mesh.tetrahedra
         for pnt in tet
            write(fid,"10 ")
         end
      end
      write(fid,"</DataArray>\n")
      write(fid,"\t\t\t</Cells>\n")

      write(fid,"\t\t\t<CellData>\n")
      for dataname in keys(data)
         write(fid,"\t\t\t\t<DataArray type=\"Float64\" Name=\"$dataname\" format=\"ascii\">")
         for datum in data[dataname]
               write(fid,string(datum)*" ")

         end
         write(fid,"\n\t\t\t\t</DataArray>\n")
      end
      write(fid,"\t\t\t</CellData>\n")
      write(fid,"\t\t</Piece>\n")
      write(fid,"\t</UnstructuredGrid>\n")
      write(fid,"</VTKFile>\n")
   end
   return nothing
end

function vtk_write2(file_name::String,mesh::Mesh,data=[])
   open(file_name*".vtu","w") do fid
      write(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" >\n")#byte_order=\"BigEndian\"
      write(fid,"\t<UnstructuredGrid>\n")
      write(fid,"\t\t<Piece NumberOfPoints=\"$(size(mesh.points,2)+length(mesh.lines))\" NumberOfCells=\"$(length(mesh.tetrahedra))\">\n")
      write(fid,"\t\t\t<Points>\n")
      write(fid,"\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")
      for pnt in mesh.points
         write(fid,string(pnt)*" ")
      end
      for line in mesh.lines
         pnts=sum(mesh.points[:,line],dims=2)/2
         for pnt in pnts
            write(fid,string(pnt)*" ")
         end
      end
      write(fid,"</DataArray>")
      write(fid,"\t\t\t</Points>\n")
      write(fid,"\t\t\t<Cells>\n")
      write(fid,"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">")
      Npoints=size(mesh.points,2)
      for tet in mesh.tetrahedra
         for pnt in tet
            write(fid,string(Int32(pnt-1))*" ")
         end
         pnt=get_line_idx(mesh,[tet[1]; tet[2]])+Npoints
         write(fid,string(Int32(pnt-1))*" ")
         pnt=get_line_idx(mesh,[tet[2]; tet[3]])+Npoints
         write(fid,string(Int32(pnt-1))*" ")
         pnt=get_line_idx(mesh,[tet[3]; tet[1]])+Npoints
         write(fid,string(Int32(pnt-1))*" ")
         pnt=get_line_idx(mesh,[tet[1]; tet[4]])+Npoints
         write(fid,string(Int32(pnt-1))*" ")
         pnt=get_line_idx(mesh,[tet[2]; tet[4]])+Npoints
         write(fid,string(Int32(pnt-1))*" ")
         pnt=get_line_idx(mesh,[tet[3]; tet[4]])+Npoints
         write(fid,string(Int32(pnt-1))*" ")
      end
      write(fid,"</DataArray>\n")
      write(fid,"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">")
      off=0
      for tet in mesh.tetrahedra
            off+=10#length(tet)
            write(fid,string(off)*" ")
      end
      write(fid,"</DataArray>\n")
      write(fid,"\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">")
      for tet in mesh.tetrahedra
         #for pnt in tet
            write(fid,"24 ")
         #end
      end
      write(fid,"</DataArray>\n")
      write(fid,"\t\t\t</Cells>\n")

      write(fid,"\t\t\t<PointData>\n")
      for dataname in keys(data)
         write(fid,"\t\t\t\t<DataArray type=\"Float64\" Name=\"$dataname\" format=\"ascii\">")
         for datum in data[dataname]
               write(fid,string(datum)*" ")

         end
         write(fid,"\n\t\t\t\t</DataArray>\n")
      end
      write(fid,"\t\t\t</PointData>\n")
      write(fid,"\t\t</Piece>\n")
      write(fid,"\t</UnstructuredGrid>\n")
      write(fid,"</VTKFile>\n")
   end
   return nothing
end
##
function vtk_write_tri(file_name::String,mesh::Mesh,data=[])
   open(file_name*".vtu","w") do fid
      write(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" >\n")#byte_order=\"BigEndian\"
      write(fid,"\t<UnstructuredGrid>\n")
      write(fid,"\t\t<Piece NumberOfPoints=\"$(size(mesh.points,2))\" NumberOfCells=\"$(length(mesh.triangles))\">\n")
      write(fid,"\t\t\t<Points>\n")
      write(fid,"\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")
      for pnt in mesh.points
         write(fid,string(pnt)*" ")
      end
      write(fid,"</DataArray>")
      write(fid,"\t\t\t</Points>\n")
      write(fid,"\t\t\t<Cells>\n")
      write(fid,"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">")
      for tri in mesh.triangles
         for pnt in tri
            write(fid,string(Int32(pnt-1))*" ")
         end
      end
      write(fid,"</DataArray>\n")
      write(fid,"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">")
      off=0
      for tri in mesh.triangles
            off+=length(tri)
            write(fid,string(off)*" ")
      end
      write(fid,"</DataArray>\n")
      write(fid,"\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">")
      for tri in mesh.triangles
         for pnt in tri
            write(fid,"5 ")
         end
      end
      write(fid,"</DataArray>\n")
      write(fid,"\t\t\t</Cells>\n")

      write(fid,"\t\t\t<CellData>\n")
      for dataname in keys(data)
         write(fid,"\t\t\t\t<DataArray type=\"Float64\" Name=\"$dataname\" format=\"ascii\">")
         for datum in data[dataname]
               write(fid,string(datum)*" ")

         end
         write(fid,"\n\t\t\t\t</DataArray>\n")
      end
      write(fid,"\t\t\t</CellData>\n")
      write(fid,"\t\t</Piece>\n")
      write(fid,"\t</UnstructuredGrid>\n")
      write(fid,"</VTKFile>\n")
   end
   return nothing
end


##

"""
    vtk_write(file_name, mesh, data)

Write vtk-file containing the datavalues `data` given on the usntructured grid `mesh`.

# Arguments
- `file_name::String`: name given to the written files. The name will be preceeded by a specifier. See Notes below.
- `mesh::Mesh`: mesh associated with the data.
- `data::Dict`: Dictionairy containing the data.

# Notes
The routine automatically sorts the data according to its type and writes it in up to three diffrent files.
Data that is constant on a tetrahedron goes into `"\$(filename)_const.vtu"`,
data that is linearly interpolated on a tetrahedron goes into `"\$(filename)_lin.vtu"`,
and data that is quadratically interpolated on a tetrahedron goes into `"\$(filename)_quad.vtu"`
"""
function vtk_write(file_name::String,mesh::Mesh,data::Dict) #TODO: write into one file
   data0=Dict()
   data1=Dict()
   data2=Dict()
   for (key,val) in data
      if length(val)==length(mesh.tetrahedra)
         data0[key]=val
      elseif length(val)==size(mesh.points,2) #TODO: check for rare case of same number of points and tets
         data1[key]=val
      elseif length(val)==size(mesh.points,2)+length(mesh.lines)
         data2[key]=val
      else
         println("WARNING: data length of '$key' does not fit to any mesh dimension.")
      end
   end

   if length(data0)!=0
      vtk_write0(file_name*"_const",mesh,data0)
   end
   if length(data1)!=0
      vtk_write1(file_name*"_lin",mesh,data1)
   end
   if length(data2)!=0
      vtk_write2(file_name*"_quad",mesh,data2)
   end
   return
end
## bloch utilities
function bloch_expand(mesh::Mesh,sol,b=:b)
    naxis=mesh.dos.naxis
    nxsector=mesh.dos.nxsector
    DOS=mesh.dos.DOS
    v=zeros(ComplexF64,naxis+nxsector*DOS)
    v[1:naxis]=sol.v[1:naxis]
    B=sol.params[b]
    for s=0:DOS-1
        v[naxis+1+s*nxsector:naxis+(s+1)*nxsector]=sol.v[naxis+1:naxis+nxsector].*exp(+2.0im*pi/DOS*B*s)
    end
    return v
end
