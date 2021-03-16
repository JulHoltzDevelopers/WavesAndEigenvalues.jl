"""
    points, lines, triangles, tetrahedra, domains=read_nastran(filename)

Read simplicial nastran mesh.

# Note
There is currently no sanity check for non-simplicial elements. These elements are just skipped while reading.
"""
function read_nastran(filename)
  number_of_points=0
  number_of_lines=0
  number_of_triangles=0
  number_of_tetrahedra=0
  name_tags=Dict()

  #first read to get number of points
  open(filename) do fid
    while !eof(fid)
      line=readline(fid)
      if line[1]=='$' #line is a comment
        #comments are freqeuntly used to leave more information
        #especially the domain tags
        #this is not standardized so it depends on the software used to write the file

        #ANSA-style tag comment
        #example:
        #$ANSA_NAME_COMMENT;2;PSHELL;TA_IN_Inlet;;NO;NO;NO;NO;
        if length(line)>=18 && line[2:18]=="ANSA_NAME_COMMENT"
          data=split(line[2:end],";")
          if data[3]=="PSOLID" || data[3]=="PSHELL"
            name_tags[data[2]]=data[4]
          end
        end

        #Centaur-style tag comment
        #example:
        #$HMNAME COMP                   1"Walls"
        if length(line)>=12 && line[2:12]=="HMNAME COMP"
          data=split(strip(line[13:end]),'"')
          name_tags[data[1]]=data[2]
        end

        continue
      end
      if length(line)<8 #line is too short probably it is the ENDDATA tag.
        continue
      end
      head=line[1:8]
      if head =="GRID    " || head == "GRID*   " || head[1:5] == "GRID," #TODO allow for whitespace between commas in free format
        number_of_points+=1
      elseif (head[1:6] == "CTRIA3") || (head[1:6] == "CTRIA6")
        number_of_triangles+=1
      elseif head[1:6] == "CTETRA"
        number_of_tetrahedra+=1
      end
    end
  end
  tri_idx=0
  tet_idx=0
  #TODO: sanity check for non-consecutive indexing
  points=Array{Float64}(undef,3,number_of_points)
  lines=Array{Array{UInt32,1}}(undef,number_of_lines)
  triangles=Array{Array{UInt32,1}}(undef, number_of_triangles)
  tetrahedra=Array{Array{UInt32,1}}(undef, number_of_tetrahedra)
  domains=Dict()
  surface_tags=Dict()
  volume_tags=Dict()
  #second read to get data
  open(filename) do fid
    while !eof(fid)
      line=readline(fid)
      if line[1]=='$' #line is a comment
        continue
      end
      # data = parse_nas_line(line)
      # if length(data)==0 #line is to short probaboly it is the ENDDATA tag.
      #   continue
      # end

      if length(line)<8 #line is too short probably it is the ENDDATA tag.
        continue
      end

      free=false #sentinel value for checking for free format
      if occursin(",",line)
        free=true
      end

      head=line[1:8]

      if head =="GRID    " #grid point short format
        data = parse_nas_line(line)
        idx=parse(UInt,data[2])
        points[1,idx]=parse_nas_number(data[4])
        points[2,idx]=parse_nas_number(data[5])
        points[3,idx]=parse_nas_number(data[6])

      elseif head =="GRID*   " #grid point long format
        data = parse_nas_line(line,format=:long)
        idx=parse(UInt,data[2])
        points[1,idx]=parse_nas_number(data[4])
        points[2,idx]=parse_nas_number(data[5])
        line=readline(fid) #long format has two lines
        data = parse_nas_line(line,format=:long)
        points[3,idx]=parse_nas_number(data[2])

      elseif head[1:5] =="GRID,"#grid point free format
        data = parse_nas_line(line,format=:free)
        idx=parse(UInt,data[2])
        points[1,idx]=parse_nas_number(data[4])
        points[2,idx]=parse_nas_number(data[5])
        points[3,idx]=parse_nas_number(data[6])

      elseif (head[1:6] == "CTRIA3") || (head[1:6] == "CTRIA6") #Triangle
        tri_idx+=1
        if free
          data = parse_nas_line(line,format=:free)
        else
          data = parse_nas_line(line)
        end
          #TODO: find nastran file specification and see whether there is
          #something like long format for stuff other than GRID
        idx=tri_idx
        dom=strip(data[3])
        if dom in keys(name_tags)
          dom=name_tags[dom]
        else
          dom="surf"*lpad(dom,4,"0")
        end
        if dom in keys(domains)
          append!(domains[dom]["simplices"],idx)
        else
          domains[dom]=Dict("dimension"=>2,"simplices"=>UInt32[idx])
        end
        tri=[parse(UInt32,data[4]),parse(UInt32,data[5]),parse(UInt32,data[6])]
        triangles[idx]=tri

      elseif head[1:6] == "CTETRA" #Tetrahedron
        tet_idx+=1
        if free
          data = parse_nas_line(line,format=:free)
        else
          data = parse_nas_line(line)
        end
        idx=tet_idx#parse(UInt,data[2])-number_of_triangles
        dom=strip(data[3])
        if dom in keys(name_tags)
          dom=name_tags[dom]
        else
          dom="vol"*lpad(dom,4,"0")
        end
        if dom in keys(domains)
          append!(domains[dom]["simplices"],idx)
        else
          domains[dom]=Dict("dimension"=>3,"simplices"=>UInt32[idx])
        end
        tet=[parse(UInt32,data[4]),parse(UInt32,data[5]),parse(UInt32,data[6]),parse(UInt32,data[7])]
        tetrahedra[idx]=tet
      end
    end
  end

  #sanity check when using higher order simplices (CTRIA6 or CTETRA with ten nodes)
  #Comsol uses this stuff -.-
  used_idx=sort(unique(vcat(tetrahedra...)))
  if used_idx[end]==length(used_idx)
    if length(triangles)!=0
      tri_max_idx=maximum(vcat(triangles...))
    else
      tri_max_idx=0
    end
    if tri_max_idx > used_idx[end]
      println("Warning: Nastran mesh needs manual reordering of point indeces. Not all points are actually used to define simplices.")
    else
      points=points[:,1:used_idx[end]]
    end
  else
    println("Warning: Nastran mesh needs manual reordering of point indeces. Not all points are actually used to define simplices.")
  end
  return points, lines, triangles, tetrahedra, domains
end


function parse_nas_line(txt;format=:short)
  #split line
  if format==:short
    number_of_entries=length(txt)÷8
    data=Array{String}(undef,number_of_entries)
    for idx=1:number_of_entries
      data[idx]=txt[(1+8*(idx-1)):(idx*8)]
    end
  elseif format==:long
    number_of_entries=(length(txt)-8)÷16+1
    data=Array{String}(undef,number_of_entries)
    data[1]=txt[1:8]
    for idx=2:number_of_entries
      data[idx]=txt[(1+16*(idx-1)-8):((idx*16)-8)]
    end
  elseif format==:free
    data=split(replace(txt," "=>""),",")
  end
  return data
end

function parse_nas_number(txt)
  # nastran numbers may or may not have sign symbols
  # and leading 0
  # also there might be lower or upper Case indicating
  # scientific notation or just a directly a sign symbol
  # This means, everything that is correctly parsed as
  # as a julia Float is a nastran float
  # but also something exotic like "-.314+7" is valid
  # nastran and would be the same as "-0.314E+7"
  # this routine correctly parses nastran numbers

  txt=strip(txt)
  E=false
  ins=0
  if !occursin("E",txt) && !occursin("e",txt)
    for (idx,c) in enumerate(txt)
      if c=="E"
        E=true
      end
      if idx!=1 && ( (!E && c=='-') || (!E && c=='+'))
        ins=idx
        break
      end
    end
    if ins >1
      txt=txt[1:(ins-1)]*"E"*txt[ins:end]
    end
  end
  return parse(Float64,txt)
end

# function relabel_nas_domains!(domains,filename)
#   open(filename) do fid
#     while !eof(fid)
#       line=readline(fid)
#       if line[1]=='#' #line is a comment
#         continue
#       end
#       data=split(line,",")
#       old_idx=strip(data[1])
#       new_idx=strip(data[2])
#       if old_idx in keys(domains)
#         domains[new_idx]=domains[old_idx]
#         delete!(domains,old_idx)
#       else
#         println("Warning: key `$old_idx` not in domains!")
#       end
#     end
#   end
#   return domains
# end
