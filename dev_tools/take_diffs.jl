
points=[(1,0,0),(0,1,0),(0,0,1),(0,0,0),(1//2,1//2,0),(1//2,0,1//2),(1//2,0,0),(0,1//2,1//2),(0,1//2,0),(0,0,1//2)]
txt="["
for (X,Y,Z) in points[1:4]
    global txt
    for func in F[2]
        df=SymPy.diff(func,z)
        num=df.subs(x,X).subs(y,Y).subs(z,Z)
        num=parse(Float64,"$num")
        txt*="$num "
    end
     txt*=";\n"
end
txt*="]"
println(txt)
