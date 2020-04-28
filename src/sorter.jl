
"""
    cmpr=compare(A,B)

Compare two simplices following lexicographic rules.
If B has lower rank, then A `cmpr ==-1`. If B has higher rank than A,
`cmpr==1`. If A and B are identical, `cmpr==0`.
"""
function compare(A,B)
    A=sort(A,rev=true)
    B=sort(B,rev=true)
    if length(B)!=length(A)
        println("Error: lists do not share same length!")
        println(A,B)
        return
    end

    cmpr=0
    for (a,b) in zip(A,B)#zip(A[end:-1:1],B[end:-1:1])
        if a!=b
            if a<b
                cmpr=1
                break
            else
                cmpr=-1
                break
            end
        end
    end
    return cmpr
end
function compare(a::Real,b::Real)
    cmpr=0
    if a!=b
        if a<b
            cmpr=1
        else
            cmpr=-1
        end
    end
    return cmpr
end



"""
    sort_smplx(list,smplx)

Helper function for sorting simplices in lexicographic manner.
Utilize divide an conquer strategy to find a simplex in ordered list of simpleces.
"""
function sort_smplx(list,smplx)
    idx=0
    list_length=length(list)
    if list_length>=2

        low=1
        high=list_length

        #check whether smplx is in list range
        # lower bound
        cmpr=compare(list[low],smplx)
        if cmpr<0
            idx=low
            return idx,true
        elseif cmpr==0
            idx=low
            return idx,false
        end

        #upper bound
        cmpr=compare(smplx,list[high])
        if cmpr<0
            idx=high+1
            return idx,true
        elseif cmpr==0
            idx=high
            return idx,false
        end

        if idx==0
            #divide et impera
            while high-low>1
                mid=Int(round((high+low)/2))

                cmpr=compare(list[mid],smplx)
                if cmpr>0
                    low=mid
                elseif cmpr<0
                    high=mid
                else
                    idx=mid
                    return idx, false
                end
            end
        end


        if idx==0
            #find location from last two possibilities
            cmpr_low=compare(list[low],smplx)
            cmpr_high=compare(smplx,list[high])
            if cmpr_low==0
                idx=low
                return idx,false
            elseif cmpr_high==0
                idx=high
                return idx,false
            else
                idx=high
                return idx,true
            end
        end


    elseif list_length==1
        cmpr=compare(list[1],smplx)
        if cmpr<0
            idx=1
            return idx, true
        elseif cmpr>0
            idx=2
            return idx,true
        else
            idx=1
            return idx,false
        end


    else list_length==0
        return 1,true
    end
    return nothing
end

"""
    insert_smplx!(list,smplx)

Insert simplex in ordered list of simplices.
"""
function insert_smplx!(list,smplx)
    #if !isa(smplx,Real)
    #    smplx=sort(smplx)
    #end
    idx,ins=sort_smplx(list,smplx)
    if ins
        insert!(list,idx,smplx)
    end
    return
end

"""
    idx = find_smplx(list,smplx)

Find index of simplex in ordered list of simplices.
If the simplex is not contained in the list, the returned index is `nothing`
"""
function find_smplx(list,smplx)
    #if !isa(smplx,Real)
    #    smplx=sort(smplx)
    #end
    idx,ins=sort_smplx(list,smplx)
    ins=!ins
    if ins
        return idx
    else
        return nothing
    end
end
