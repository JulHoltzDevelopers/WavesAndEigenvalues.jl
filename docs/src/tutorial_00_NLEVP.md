```@meta
EditURL = "<unknown>/tutorial_00_NLEVP.jl"
```

# Tutorial 00 Introduction to the NLEVP module

This tutorial briefly introduces you to the NLEVP module.
NLEVP is an abbreviation for nonlinear eigenvalue problem.
This is the problem `T(λ)v=0` where T is some matrix family
the scalar `λ` and the vector `v` form an eigenpair. Often also the left
eigenvector `w` for the corresponding problem `T'(λ)w=0` is also sought.
Clearly, NLEVPs are a generalization of classic (linear) eigenvalue problems
they would read `T(λ)=A-λI` or `T(λ)=A-λB`. For the reminder of the tutorial
it is assumed, that the reader is somewhat familiar with NLEVPs and we refer
to the very nice review [1] for a detailed introduction.


## Setting up a NLEVP

The NLEVP lets you define NLEVPs and provides you with tools for finding
accurate as well as perturbative solutions. The `using` statement brings all
necessary tools into scope:

```julia
using WavesAndEigenvalues.NLEVP
```

Let's consider the quadratic eigenvalue problem  1 from the collection of NLEVPs [2].
It reads `T(λ)= λ^2*A2 + λ*A1 + A0` where

```julia
A2=[0 6 0;
    0 6 0;
    0 0 1];

A1=[1 -6 0;
    2 -7 0;
    0  0 0];

A0=[1 0 0;
    0 1 0;
    0 0 1];
nothing #hide
```

We first need to instantiate an empty operator family.

```julia
T=LinearOperatorFamily()
```

```
empty operator family



Parameters
----------
λ	NaN + NaN*im

```

!!! tip
    Per default the LinearOperatorFamily uses `:λ` as the symbol for
    the eigenvalue. You can customize this behavior by calling the constructor
    with additional parameters, e.g. `T=LinearOperatorFamily([:ω])`, would make
    the eigenvalue name `:ω`. It is also possible to change the eigenvalue after
    construction.

The next step is to create the three terms and add them to the matrix family.
A term is created from the syntax `Term(M,func,args,txt,mat)`, where `M` is the
matrix of the term, `func` is a list of scalar (possibly multi-variate)
functions multiplying, `args`, is a list of lists containing the symbols of the
arguments to the functions, and `mat` is the character string representing the
matrix.

Ok, let's create the first term:

```julia
term=Term(A2,(pow2,),((:λ,),),"A2")
```

```
λ^2*A2
```

Note, that the function `pow2` is provided by the NLEVP module. It computes
the square of its input, i.e. pow2(3)==9. However, it has an optional second
parameter indicating a derivative order. For instance, `pow2(3,1)=6` because
the first derivative of the square function at `3` is `2*3==6`. This feature
is, important because many algorithms in the NLEVP module will need to compute
derivatives. You can write your own coefficient-functions, but make sure
that they provide at least the first derivative correctly. This will be
further explained in the next tutorials.

Now, let's add the term to our family.

```julia
T+=term
```

```
3×3-dimensional operator family: 

λ^2*A2

Parameters
----------
λ	NaN + NaN*im

```

It worked! The family is not empty anymore.

!!! tip
    There is also the syntax `push!(T,term)` to add a term. It is
    actually more efficient because unlike the `+`-operator it does not create
    a temporary copy of `T`. However, it doesn't really matter for small problems
    like this example and we therefore prefer the `+`-operator for readability.

Let's also add the other two terms.

```julia
T+=Term(A1,(pow1,),((:λ,),),"A1")
T+=Term(A0,(),(),"A0")
```

```
3×3-dimensional operator family: 

λ^2*A2+λ*A1+A0

Parameters
----------
λ	NaN + NaN*im

```

Note that functions and arguments are empty in the last term because there is
no coefficient function.

The family is now a complete representation of our problem

It is quite powerful. For instance, we can evaluate it at some point in the
complex plane, say 3+2im:

```julia
T(3+2im)
```

```
3×3 Array{Complex{Float64},2}:
 4.0+2.0im  12.0+60.0im  0.0+0.0im
 6.0+4.0im  10.0+58.0im  0.0+0.0im
 0.0+0.0im   0.0+0.0im   6.0+12.0im
```

We can even take its derivatives. For instance, the second derivative at the
same point is

```julia
T(3+2im,2)
```

```
3×3 Array{Complex{Float64},2}:
 0.0+0.0im  12.0+0.0im  0.0+0.0im
 0.0+0.0im  12.0+0.0im  0.0+0.0im
 0.0+0.0im   0.0+0.0im  2.0+0.0im
```

This syntax comes in handy when writing algorithms for the solution of the
NLEVPs. However, the casual user will appreciate it when computing residuals.

## Iterative eigenvalue solvers

The NLEVP module provides various solvers for computing eigensolution. They
essentially fall into two categories, *iterative* and *integration-based*. We
will start the discussion with iterative solvers.

### method of successive linear problems

The method of successive linear problems is a robust iterative solver.
Like all iterative solvers it tries finding an eigensolution from an
initial guess. We will also set the `output` keyword to true, in order to see
a bit what is going on behind the scenes. The function will return a solution
`sol` of type `Solution`, as well as the number of iterations `n` run and an
error flag `flag`.

```julia
sol,n,flag=mslp(T,0,output=true);
nothing #hide
```

```
Launching MSLP solver...
Iter   dz:     z:
----------------------------------
0			Inf	0.0
1			0.21976270357743133	0.2147594636259946 + 0.04662637308151951im
2			0.08885531760332352	0.2993866118413472 + 0.019542923313732537im
3			0.033202083130972414	0.3292202626020246 + 0.004971320569835911im
4			0.006293004538713664	0.3333645192058974 + 0.0002356006324317284im
5			0.00023769944078740377	0.33333366091326566 - 8.727546025417629e-8im
6			3.390074774396457e-7	0.33333333333273524 + 3.430782475273895e-13im
7			6.893955075745975e-13	0.3333333333333332 + 2.4613974855833593e-24im
8			6.661338147750939e-16	0.33333333333333254 + 2.9321105043677363e-31im
9			6.661338147750939e-16	0.3333333333333332 + 3.339320211059557e-31im
10		6.661338147750932e-16	6.661338147750939e-16	0.33333333333333254 + 4.351492030582117e-31im
Warning: Maximum number of iterations has been reached!
...finished MSLP!
#####################
 Results 
#####################
Number of steps: 10
Last step parameter variation:6.661338147750939e-16
Auxiliary eigenvalue __aux__ residual (rhs):6.661338147750932e-16
Eigenvalue:0.33333333333333254 + 4.351492030582117e-31im

```

The printed output tells us that the algorithm found 1/3 is an eigenvalue of
the problem. It has run for 10 iterations. However, it issued a warning that
the maximum number of iterations has been reached. That is why the result
comes with a non-zero flag, which says that something might be wrong.
In general negative flag values represent errors, while positive values
indicate warnings. In the current case the flag is `1`, so this is a warning.
Therefore, the result can be true, but there is something we should be aware
of.  Let's decode the error flag into something more readable.

```julia
msg = decode_error_flag(flag)
println(msg)
```

```
Warning: Maximum number of iterations has been reached!

```

As we already knew, the iteration aborted because the maximum number of
iterations (10) has been performed. You can change the maximum number of
iterations from its default value by using the `maxiter` keyword. However,
this won't fix the problem because the iteration will run until a certain
tolerance is met or the maximum number of iteration is reached. Per default
the tolerance is 0, and therefore due to machine-precision very unlikely
to be precisely met. Indeed, we see from the printed output, above that the
eigenvalue is not significantly changing anymore after 6 iterations. Let's
redo the calculations, with a small but non-zero tolerance by setting the
optional `tol` keyword.

```julia
sol,n,flag=mslp(T,0,tol=1E-10,output=true);
nothing #hide
```

```
Launching MSLP solver...
Iter   dz:     z:
----------------------------------
0			Inf	0.0
1			0.3091753342130375	0.3084490799844025 - 0.02117905433486341im
2			0.028554443587410515	0.3316559062828601 - 0.004541763036703544im
3			0.004788776681260656	0.33343303494632615 - 9.494441925653792e-5im
4			0.00013775902065808008	0.33333332790060327 + 1.1373560610378447e-7im
5			1.1386527920077728e-7	0.33333333333341086 + 7.414843609558661e-15im
6		7.762646518664894e-14	7.762646518663691e-14	0.3333333333333336 - 6.811025055678216e-27im
Solution has converged!
...finished MSLP!
#####################
 Results 
#####################
Number of steps: 6
Last step parameter variation:7.762646518663691e-14
Auxiliary eigenvalue __aux__ residual (rhs):7.762646518664894e-14
Eigenvalue:0.3333333333333336 - 6.811025055678216e-27im

```

Boom! The iteration stops after 6 iterations and indicates that an
eigenvalue has been found.

We may also inspect the corresponding left and right eigenvector. They are
fields in the solution object:

```julia
println(sol.v)
println(sol.v_adj)
```

```
Complex{Float64}[-0.2893267018214533 + 0.6452054398508432im, -0.2893267018214534 + 0.6452054398508428im, 1.2325951644078312e-32 - 6.162975822039156e-33im]
Complex{Float64}[-0.5786534036430183 + 1.2904108797020055im, 1.1573068072859303 - 2.5808217594036975im, -1.3382229383269777e-43 + 8.938109037958203e-44im]

```

There are more iterative solvers. Not all of them solve the complete problem.
Some only return the eigenvalue and the right eigenvector others even just the
eigenvalue. (If the authors of this software get more funding this might be
improved in the future...)

Let's briefly test them all...

### inverse iteration

This algorithm finds an eigenvalue and a corresponding right eigenvector

```julia
sol,n,flag=inveriter(T, 0, tol=1E-10, output=true);
nothing #hide
```

```
Launching inverse iteration...
0		Inf	0
1		0.30000000000000004	0.30000000000000004 + 0.0im
2		0.028571428571428137	0.3285714285714282 + 0.0im
3		0.00463320463320499	0.33320463320463317 + 0.0im
4		0.00012860089961053145	0.3333332341042437 + 0.0im
5		9.922903054793153e-8	0.33333333333327425 + 0.0im
6		5.950795411990839e-14	0.33333333333333376 + 0.0im
Solution has converged!

```

### trace iteration
This algorithm only finds an eigenvalue

```julia
sol,n,flag=traceiter(T, 0, tol=1E-10, output=true);
nothing #hide
```

```
Launching trace iteration...
0		Inf	0
1		0.16666666666666666	0.16666666666666666 + 0.0im
2		0.10125889436234267	0.26792556102900933 + 0.0im
3		0.04886705474089864	0.31679261576990797 + 0.0im
4		0.014969376489607611	0.3317619922595156 + 0.0im
5		0.0015546254085949673	0.33331661766811055 + 0.0im
6		1.6713737663487382e-5	0.33333333140577404 + 0.0im
7		1.927559389880429e-9	0.3333333333333334 + 0.0im
8		2.220446049250313e-16	0.3333333333333332 + 0.0im
Solution has converged!

```

### Lancaster's generalized Rayleigh quotient iteration
This algorithm only finds an eigenvalue

```julia
sol,n,flag=lancaster(T, 0, tol=1E-10, output=true);
nothing #hide
```

```
Launching Lancaster's Rayleigh-quotient iteration...
0		Inf	0
1		0.30000000000000004	0.30000000000000004 + 0.0im
2		0.029104073704340316	0.32910407370434036 + 0.0im
3		0.004135217444984907	0.33323929114932527 + 0.0im
4		9.399316044850226e-5	0.33333328430977377 + 0.0im
5		4.902354600044845e-8	0.33333333333331977 + 0.0im
6		1.354472090042691e-14	0.3333333333333333 + 0.0im
Solution has converged!

```

### two-sided Rayleigh-functional iteration
This algorithm finds a complete eigentriple.

```julia
sol,n,flag=rf2s(T, 0, tol=1E-10, output=true);
nothing #hide
```

```
Launching two-sided Rayleigh functional iteration...
0		Inf	0
1		0.23434529628710668	0.23434529628710668 - 0.0im
2		0.09576672895627122	0.3301120252433779 - 0.0im
3		0.0032212122529545195	0.3333332374963324 - 0.0im
4		9.583700116833072e-8	0.3333333333333336 - 0.0im
5		NaN	NaN + NaN*im
Warning: computer arithmetics problem. Eigenvalue is NaN

```

Oops, this last test produced a NaN value. The problem often occurs when
the eigenvalue and one point in the iteration is too precise. (See how in the
last-but-one iteration we've been already damn close to 1/3). Choosing a
slightly different initial guess often leverages the problem:

```julia
sol,n,flag=rf2s(T, 0, tol=1E-10, output=true);
nothing #hide
```

```
Launching two-sided Rayleigh functional iteration...
0		Inf	0
1		0.23434529628710668	0.23434529628710668 - 0.0im
2		0.09576672895627122	0.3301120252433779 - 0.0im
3		0.0032212122529545195	0.3333332374963324 - 0.0im
4		9.583700116833072e-8	0.3333333333333336 - 0.0im
5		NaN	NaN + NaN*im
Warning: computer arithmetics problem. Eigenvalue is NaN

```

## Beyn's integration-based solver

Iterative solvers are highly accurate and do not need many resources. However,
they only find one eigenvalue at the time and heavily depend on the initial
guess. This is a problem because NLEVPs have many eigenvalues. (Indeed, some
have infinitely many!) This is where integration-based solvers enter the stage.
They can find all eigenvalue enclosed by a contour in the complex plane.
Several such algorithm exist. The NLEVP module implements the one of Beyn [3]
(both of its variants!). Any polygonal shape is supported as a contour
(I am an engineer, so a polygon is any flat shape connecting three or more
vertices with straight lines in a cyclic manner...).

Let's test Beyn's algorithm on our example problem by searching for all
eigenvalues inside the axis-parallel square spanning from `-2-2im` to `2+2im`.
The vertices of this contour are:

```julia
Γ=[2+2im, -2+2im, -2-2im, +2-2im];
nothing #hide
```

It doesn't matter whether the contour is traversed in clockwise or
counter-clockwise direction. However, the vertices should be traversed in a
cyclic manner such that the contour is not criss-crossed (unless this is what
you want to do, the algorithm correctly accounts for the winding number...)

Beyn's algorithm also needs an educated guess on how many eigenvalues you are
expecting to find inside of the contour. If your guess is too high, this won't
be a problem. You will find your eigenvalues but you'll waste computational
resources and face an unnecessarily high run time. Nevertheless, if your
guess is too small, weird things may happen and you probably will get wrong
results. We will further discuss this point later. Now let's focus on the
example problem. It is 3-dimensional and quadratic, we therefore know that it
that it has 6 eigenvalues. Not all of them must lie in our contour, however
this upper bound is a safe guess for our contour. We use the optional keyword
`l=6` to pass this information to the algorithm. The call then is:

```julia
Λ,V=beyn(T,Γ,l=6, output=true);
nothing #hide
```

```
Beyn...  27%|████████                      |  ETA: 0:00:02Beyn... 100%|██████████████████████████████| Time: 0:00:00
############
singular values:
[14.426760137347847, 6.283184515316616, 6.283184515316615, 3.4297851239318184, 0.6821925100199037, 1.4050439819523788e-15]

```

Note the printed information on singular values. One step in Beyn's algorithm
leads to a singular value decomposition. Here, we got 6 singular values
because we chose to run the algorithm with `l=6` option. In exact arithmetics
There would be as many non-zero singular values as there are eigenvalues
inside the contour, if the guess `l` is equal or greater than the number of
eigenvalues. The output tells us that there are probably 5 eigenvalues inside
the contour as there is only one singular value that can be considered equal
to `0` given the machine precision. Indeed, the example solution is
analytical solvable and it is known that it has 5 eigenvalues inside the
square.

The return values `Λ` and `V` are the list of eigenvalues and corresponding
(right) eigenvectors such that `T(Λ[i])*V[:,i]==0`. At least this is what
it should be. Because the computation is carried out on a computer the
eigenpairs won't evaluate to the zero vector but have some residual. We can
use this as a quality check. The computation of the residual norms yields:

```julia
for (idx,λ) in enumerate(Λ)
   v=V[:,idx]
   v/=sqrt(v'*v)
   res=T(λ)*v
   println("λ=$λ residual norm: $(abs(sqrt(res'*res)))")
end
```

```
λ=-3.0046846037705304e-16 + 1.0000000000000002im residual norm: 3.2378127103262535e-15
λ=1.9081958235744878e-16 - 0.9999999999999998im residual norm: 1.9230280355021297e-15
λ=0.16554342477066233 - 1.0271955541799844im residual norm: 3.8796425923748163
λ=0.33333333333333576 + 6.677193806922799e-16im residual norm: 2.1182915745378617e-15
λ=0.49999999999999545 - 1.0637329323545343e-15im residual norm: 2.528153131782494e-15
λ=1.000000000000001 - 1.5080235779736662e-16im residual norm: 6.476979535860806e-15

```

We clearly see that out of the 6 eigenvalues 5 have extremely low residuals.
The eigenvalue with the high residual is actually not an eigenvalue. This is
inline with our consideration on the singular value. There is a keyword
`sigma_tol` that you may provide to the call of the algorithm, in order to
ignore singular values that are smaller then a user-specified threshold.# However, an appropriate threshold is problem-dependent and often it is a
better strategy to feed the computed eigenvalues as initial guesses into an
iterative solver.

Please, note that you could also compute the number of poles and zeros inside the contour from
the residual theorem with the following command:

```julia
number=count_poles_and_zeros(T,Γ)
```

```
4.999999621784568 - 3.092205940311753e-17im
```

The functions sums up the number of poles and zeros of the determinant of `T`
inside your contour. Zeros are positively counted while poles count negatively
when encircling the domain mathematically positive (counter-clockwise). Thus,
the number should not be mistaken as an initial guess for `l`. Yet, often
you have some knowledge about the properties of your problem. For instance
the current example problem is polynomial an therefore is unlikely to have
poles. That `number≈5` is a further indicator that there are, indeed, only 5
eigenvalues inside the square. Note, that `count_poles_and_zeros` uses Jacobi's
formula to compute the necessary derivatives of the determinant by trace
operations. The command is therefore reliable for small-sized problems only.


## Summary

The NLEVP module provides you with essential tools for setting up and solving
non-linear eigenvalue problems. In particular, there is a bunch of solution
algorithms available. There is more the module can do. All of the presented
solvers have various optional keyword arguments for fine-tuning there
behavior. For instance, just type `? mslp` in the REPL or use the
documentation browser of your choice to learn more about the `mslp` command.
Moreover, there is functionality for computing perturbation expansions of
known solutions. These topics will be further discussed in the next tutorials.

# References

[1] S. Güttel and F. Tisseur, The Nonlinear Eigenvalue Problem, 2017, <http://eprints.ma.man.ac.uk/2538/>

[2] T. Betcke, N. J. Higham, V. Mehrmann, C. Schröder and F. Tisseur, NLEVP: A Collection of Nonlinear Eigenvalue Problems, 2013, ACM Trans. Math. Soft., [doi:10.1145/2427023.2427024](https://doi.org/10.1145/2427023.2427024)

[3] W.-J. Beyn, An integral method for solving nonlinear eigenvalue problems, 2012, Lin. Alg. Appl., [doi:10.1016/j.laa.2011.03.030](https://doi.org/10.1016/j.laa.2011.03.030)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

