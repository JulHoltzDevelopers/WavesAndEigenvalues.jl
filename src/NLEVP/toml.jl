#module TOML
#export read_toml
#TODO: tags cant be identical to variable names
#TODO: optimize parsing of multiline lists
#TODO: nested data
#TODO: preallocation?
#TODO: put this in a module to correctly captur __data_TOML__
#TODO: write a TOPL type to conveniently access the dictionary with TOML-like syntax
"This is a minimal TOML parser. See https://en.wikipedia.org/wiki/TOML"
function read_toml(fname)
    tags=""
    multi=false #sentinel value for multiline input
    data,var="","" #make data and var global in loop
    D=Dict()
    entry=Dict() # make entry global in loop
    for (linenumber,line) in enumerate(eachline(fname))
        line=strip(line)#remove enclosing whitespace
        if isempty(line) || line[1]=="#"
            continue
        elseif !multi && line[1]=='['
            #tag
            #TODO: check whether last chracter is ] otherwise error
            #tags=split(strip(line,['[',']']),'.') #parse tag
            #TODO: remove enclosing intermediate whitespace?
            tags=split(line[2:end-1],'.')
            entry=D
            for tag in tags
                tag="/"*tag
                if tag in keys(entry)
                    entry=entry[tag]
                else
                    #TODO: consider signaling an error if tag is not the last entry in tags
                    entry[tag]=Dict()
                    entry=entry[tag]
                end
            end
        elseif !multi && isletter(line[1])
            idx=findfirst("=",line)
            var=line[1:idx.start-1]
            var=strip(var) #remove enclosing whitespace
            data=line[idx.stop+1:end]
            multi=data[end]==','
        elseif multi
            data*=line
            multi=data[end]==','
        else
            #TODO:Throw error
        end

        if !multi && data!=""
            #println(data)
            eval(Meta.parse("__data_TOML__="*data)) #evaluate data in global scope of module
            if tags ==""
                D[var]=__data_TOML__
            else
                entry[var]=__data_TOML__
            end
            data=""
        end
    end
    eval(Meta.parse("__data_TOML__=nothing"))
    return D
end
#end
