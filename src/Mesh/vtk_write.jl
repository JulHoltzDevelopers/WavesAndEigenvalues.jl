#functionality in this file implements writing .vtu files
#complient to https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
import Base64
#import CodecZlib
##
function vtk_write1(file_name::String,mesh::Mesh,data=[],binary=true,compressed=false)
   open(file_name*".vtu","w") do fid
      if compressed
         compressor=" compressor=\"vtkZLibDataCompressor\""
      else
         compressor=""
      end
      write(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\"$compressor>\n")
      write(fid,"\t<UnstructuredGrid>\n")
      write(fid,"\t\t<Piece NumberOfPoints=\"$(size(mesh.points,2))\" NumberOfCells=\"$(length(mesh.tetrahedra))\">\n")
      write(fid,"\t\t\t<Points>\n")
      #write(fid,"\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")
      #for pnt in mesh.points
      #   write(fid,string(pnt)*" ")
      #end
      #write(fid,"</DataArray>")
      write_data_array(fid,Dict("Coordinates"=>mesh.points),binary,compressed)
      write(fid,"\t\t\t</Points>\n")
      write(fid,"\t\t\t<Cells>\n")
      # write(fid,"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">")
      # for tet in mesh.tetrahedra
      #    for pnt in tet
      #       write(fid,string(Int32(pnt-1))*" ")
      #    end
      # end
      # write(fid,"</DataArray>\n")
      N_points=size(mesh.points,2)
      if N_points <2^16
         connectivity=Array{UInt16}(undef,length(mesh.tetrahedra)*4)
         idx=1
         for tet in mesh.tetrahedra
            for pnt_idx in tet
               connectivity[idx]=UInt16(pnt_idx-1)
               idx+=1
            end
         end
      else
         connectivity=Array{UInt32}(undef,length(mesh.tetrahedra)*4)
         idx=1
         for tet in mesh.tetrahedra
            for pnt_idx in tet
               connectivity[idx]=UInt32(pnt_idx-1)
               idx+=1
            end
         end
      end
      write_data_array(fid,Dict("connectivity"=>connectivity),binary,compressed)

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
      write_data_array(fid,data,binary,compressed)
      write(fid,"\t\t\t</PointData>\n")
      write(fid,"\t\t</Piece>\n")
      write(fid,"\t</UnstructuredGrid>\n")
      write(fid,"</VTKFile>\n")
   end
   return nothing
end
##
function vtk_write0(file_name::String,mesh::Mesh,data=[],binary=true,compressed=false)
   open(file_name*".vtu","w") do fid
      if compressed
         compressor=" compressor=\"vtkZLibDataCompressor\""
      else
         compressor=""
      end
      write(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\"$compressor>\n")
      write(fid,"\t<UnstructuredGrid>\n")
      write(fid,"\t\t<Piece NumberOfPoints=\"$(size(mesh.points,2))\" NumberOfCells=\"$(length(mesh.tetrahedra))\">\n")
      write(fid,"\t\t\t<Points>\n")
      # write(fid,"\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")
      # for pnt in mesh.points
      #    write(fid,string(pnt)*" ")
      # end
      # write(fid,"</DataArray>")
      write_data_array(fid,Dict("Coordinates"=>mesh.points),binary,compressed)
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
      write_data_array(fid,data,binary,compressed)
      write(fid,"\t\t\t</CellData>\n")
      write(fid,"\t\t</Piece>\n")
      write(fid,"\t</UnstructuredGrid>\n")
      write(fid,"</VTKFile>\n")
   end
   return nothing
end

function vtk_write2(file_name::String,mesh::Mesh,data=[],binary=true,compressed=false)
   open(file_name*".vtu","w") do fid
      if compressed
         compressor=" compressor=\"vtkZLibDataCompressor\""
      else
         compressor=""
      end
      write(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\"$compressor>\n")
      write(fid,"\t<UnstructuredGrid>\n")
      write(fid,"\t\t<Piece NumberOfPoints=\"$(size(mesh.points,2)+length(mesh.lines))\" NumberOfCells=\"$(length(mesh.tetrahedra))\">\n")
      write(fid,"\t\t\t<Points>\n")
      #write(fid,"\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")
      # for pnt in mesh.points
      #    write(fid,string(pnt)*" ")
      # end
      # for line in mesh.lines
      #    pnts=sum(mesh.points[:,line],dims=2)/2
      #    for pnt in pnts
      #       write(fid,string(pnt)*" ")
      #    end
      # end
      #write(fid,"</DataArray>")
      N_points=size(mesh.points,2)
      points=Array{Float64}(undef,3,N_points+length(mesh.lines))
      points[:,1:N_points]=mesh.points
      for (idx,line) in enumerate(mesh.lines)
         points[:,N_points+idx]=sum(mesh.points[:,line],dims=2)/2
      end
      write_data_array(fid,Dict("coordinates"=>points),binary,compressed)
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
      write_data_array(fid,data,binary,compressed)
      write(fid,"\t\t\t</PointData>\n")
      write(fid,"\t\t</Piece>\n")
      write(fid,"\t</UnstructuredGrid>\n")
      write(fid,"</VTKFile>\n")
   end
   return nothing
end
##
function vtk_write_tri(file_name::String,mesh::Mesh,data=[],binary=true,compressed=false)
   open(file_name*".vtu","w") do fid
      if compressed
         compressor=" compressor=\"vtkZLibDataCompressor\""
      else
         compressor=""
      end
      write(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\"$compressor>\n")
      write(fid,"\t<UnstructuredGrid>\n")
      write(fid,"\t\t<Piece NumberOfPoints=\"$(size(mesh.points,2))\" NumberOfCells=\"$(length(mesh.triangles))\">\n")
      write(fid,"\t\t\t<Points>\n")
      # write(fid,"\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")
      # for pnt in mesh.points
      #    write(fid,string(pnt)*" ")
      # end
      # write(fid,"</DataArray>")
      write_data_array(fid,Dict("Coordinates"=>mesh.points),binary,compressed)
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
      write_data_array(fid,data,binary,compressed)
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
- `data::Dict`: Dictionary containing the data.

# Notes
The routine automatically sorts the data according to its type and writes it in up to three diffrent files.
Data that is constant on a tetrahedron goes into `"\$(filename)_const.vtu"`,
data that is linearly interpolated on a tetrahedron goes into `"\$(filename)_lin.vtu"`,
and data that is quadratically interpolated on a tetrahedron goes into `"\$(filename)_quad.vtu"`
"""
function vtk_write(file_name::String,mesh::Mesh,data::Dict,binary=false,compressed=false) #TODO: write into one file
   data0=Dict()
   data1=Dict()
   data2=Dict()
   data_tri=Dict()
   for (key,val) in data
      if maximum(size(val))==length(mesh.tetrahedra)
         data0[key]=val
      elseif maximum(size(val))==size(mesh.points,2) #TODO: check for rare case of same number of points and tets
         data1[key]=val
      elseif maximum(size(val))==size(mesh.points,2)+length(mesh.lines)
         data2[key]=val
      elseif maximum(size(val))==length(mesh.triangles)
         data_tri[key]=val
      else
         println("WARNING: data length of '$key' does not fit to any mesh dimension.")
      end
   end

   if length(data0)!=0
      vtk_write0(file_name*"_const",mesh,data0,binary,compressed)
   end
   if length(data1)!=0
      vtk_write1(file_name*"_lin",mesh,data1,binary,compressed)
   end
   if length(data2)!=0
      vtk_write2(file_name*"_quad",mesh,data2,binary,compressed)
   end
   if length(data_tri)!=0
      vtk_write_tri(file_name*"_tri",mesh,data_tri,binary,compressed)
   end
   return
end
## bloch utilities
# function bloch_expand(mesh::Mesh,sol,b=:b)
#     naxis=mesh.dos.naxis
#     nxsector=mesh.dos.nxsector
#     DOS=mesh.dos.DOS
#     v=zeros(ComplexF64,naxis+nxsector*DOS)
#     v[1:naxis]=sol.v[1:naxis]
#     B=sol.params[b]
#     for s=0:DOS-1
#         v[naxis+1+s*nxsector:naxis+(s+1)*nxsector]=sol.v[naxis+1:naxis+nxsector].*exp(+2.0im*pi/DOS*B*s)
#     end
#     return v
# end
# function bloch_expand(mesh::Mesh,vec::Array,B::Real=0)
#     naxis=mesh.dos.naxis
#     nxsector=mesh.dos.nxsector
#     DOS=mesh.dos.DOS
#     v=zeros(ComplexF64,naxis+nxsector*DOS)
#     v[1:naxis]=sol.v[1:naxis]
#     for s=0:DOS-1
#         v[naxis+1+s*nxsector:naxis+(s+1)*nxsector]=vec[naxis+1:naxis+nxsector].*exp(+2.0im*pi/DOS*B*s)
#     end
#     return v
# end

## helper functions for writing vtk data

function write_data_array(fid,data,binary=true,compressed=false)
   if binary
      if compressed
         iob = IOBuffer()
         #io.append=true
         #iob = Base64.Base64EncodePipe(iob);
      else
         iob = Base64.Base64EncodePipe(fid);
      end
      format="\"binary\""
   else
      format="\"ascii\""
   end

   for dataname in keys(data)
      if size(data[dataname],2)>1 #write vector output #TODO: Sanity checks for three components
         datatype="$(typeof(data[dataname][1]))"
         bytesize=sizeof(data[dataname][1])
         write(fid,"\t\t\t\t<DataArray type=\"$datatype\" NumberOfComponents=\"3\" Name=\"$dataname\" format=$format>")
         write(fid,"\n\t\t\t\t\t")
         if binary
            len=length(data[dataname])
            if !compressed
               num_of_bytes=UInt64(length(data[dataname])*bytesize)
               write(iob, num_of_bytes)
            end
            for vec in data[dataname]
               for num in vec
                  write(iob, num)#write number to base64 coded stream
               end
            end
         else
            for vec in data[dataname]
               write(fid,string(vec)*" ")
            end
         end
      else
         datatype="$(typeof(data[dataname][1]))"
         bytesize=sizeof(data[dataname][1])
         #datatype="Float64"
         write(fid,"\t\t\t\t<DataArray type=\"$datatype\" Name=\"$dataname\" format=$format>")
         write(fid,"\n\t\t\t\t\t")
         if binary
               #num_of_bytes=UInt64(length(data[dataname])*bytesize)
            if !compressed
               num_of_bytes=UInt64(length(data[dataname])*bytesize)
               write(iob, num_of_bytes)
            end
            for datum in data[dataname]
               idxx=0
               #println("$idxx : $dataname")
               for (idxc,num) in enumerate(datum)
                  idxx+=1
                  #println("$idxx : $num")
                  write(iob, num)#write number to base64 coded stream
               end
            end
         else
            for datum in data[dataname]
                  write(fid,string(datum)*" ")
            end
         end
      end
      if binary
         #close(iob) #close io stream
         #zip compress data
         if compressed
            block_zcompress!(fid,take!(iob))
         end
         #write(fid,String(take!(io))) #write encoded data to file
      end
      write(fid,"\n\t\t\t\t</DataArray>\n")
   end
   if compressed
      close(iob)
   end
end

function block_zcompress!(fid,data,block_length=2^10)
    iob64_encode = Base64.Base64EncodePipe(fid);
    len=length(data)
    number_of_blocks=len÷block_length
    length_of_last_partial_block=len%block_length
    if length_of_last_partial_block!=0
        number_of_blocks+=1
    end
    compressed_block_size = Array{UInt64}(undef,number_of_blocks)
    compressed_data=[]
    for idx=1:number_of_blocks
        if idx==number_of_blocks
            block=data[((idx-1)*block_length+1):end]
        else
            block=data[((idx-1)*block_length+1):(idx*block_length)]
        end
        compressed_block = CodecZlib.transcode(CodecZlib.GzipCompressor, block)
        compressed_block_size[idx]=sizeof(compressed_block)
        append!(compressed_data,[compressed_block])
    end

    #write compressed binary entry
    #io = IOBuffer()
    #header
    write(iob64_encode, UInt64(number_of_blocks))
    write(iob64_encode, UInt64(block_length))
    write(iob64_encode, UInt64(length_of_last_partial_block))
    for idx=1:number_of_blocks
      write(iob64_encode, UInt64(compressed_block_size[idx]))
    end
    #data
    for idx=1:number_of_blocks
        for num in compressed_data[idx]
            write(iob64_encode, num)
        end
    end

    return nothing
end
