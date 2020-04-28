#include("mesh.jl")
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
end
##
#data=Dict()
#data[string(round(ω/2/pi,2))]=abs.(p)
#vtk_write("mytest",mesh,data)


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
         pnt=find_smplx(mesh.lines,[tet[1]; tet[2]])+Npoints
         write(fid,string(Int32(pnt-1))*" ")
         pnt=find_smplx(mesh.lines,[tet[2]; tet[3]])+Npoints
         write(fid,string(Int32(pnt-1))*" ")
         pnt=find_smplx(mesh.lines,[tet[3]; tet[1]])+Npoints
         write(fid,string(Int32(pnt-1))*" ")
         pnt=find_smplx(mesh.lines,[tet[1]; tet[4]])+Npoints
         write(fid,string(Int32(pnt-1))*" ")
         pnt=find_smplx(mesh.lines,[tet[2]; tet[4]])+Npoints
         write(fid,string(Int32(pnt-1))*" ")
         pnt=find_smplx(mesh.lines,[tet[3]; tet[4]])+Npoints
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
end

function vtk_write(file_name::String,mesh::Mesh,data=[]) #TODO: write into one file
   datum = first(values(data))
   if length(datum)==length(mesh.tetrahedra)
      vtk_write0(file_name,mesh,data)
   elseif length(datum)==size(mesh.points,2) #TODO: check for rare case of same number of points and tets
      vtk_write1(file_name,mesh,data)
   elseif length(datum)==size(mesh.points,2)+length(mesh.lines)
      vtk_write2(file_name,mesh,data)
   else
      println("ERROR: data length does not fit to any mesh dimension")
   end
   return
end
