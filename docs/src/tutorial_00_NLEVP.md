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

```@example tutorial_00_NLEVP
using WavesAndEigenvalues.NLEVP
```

Let's consider the quadratic eigenvalue problem  1 from the collection of NLEVPs [2].
It reads `T(λ)= λ^2*A2 + λ*A1 + A0` where

```@example tutorial_00_NLEVP
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

```@example tutorial_00_NLEVP
T=LinearOperatorFamily()
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

```@example tutorial_00_NLEVP
term=Term(A2,(pow2,),((:λ,),),"A2")
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

```@example tutorial_00_NLEVP
T+=term
```

It worked! The family is not empty anymore.

!!! tip
    There is also the syntax `push!(T,term)` to add a term. It is
    actually more efficient because unlike the `+`-operator it does not create
    a temporary copy of `T`. However, it doesn't really matter for small problems
    like this example and we therefore prefer the `+`-operator for readability.

Let's also add the other two terms.

```@example tutorial_00_NLEVP
T+=Term(A1,(pow1,),((:λ,),),"A1")
T+=Term(A0,(),(),"A0")
```

Note that functions and arguments are empty in the last term because there is
no coefficient function.

The family is now a complete representation of our problem

It is quite powerful. For instance, we can evaluate it at some point in the
complex plane, say 3+2im:

```@example tutorial_00_NLEVP
T(3+2im)
```

We can even take its derivatives. For instance, the second derivative at the
same point is

```@example tutorial_00_NLEVP
T(3+2im,2)
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

```@example tutorial_00_NLEVP
sol,n,flag=mslp(T,0,output=true);
nothing #hide
```

The printed output tells us that the algorithm found 1/3 is an eigenvalue of
the problem. It has run for 10 iterations. However, it issued a warning that
the maximum number of iterations has been reached. That is why the result
comes with a non-zero flag, which says that something might be wrong.
In general negative flag values represent errors, while positive values
indicate warnings. In the current case the flag is `1`, so this is a warning.
Therefore, the result can be true, but there is something we should be aware
of.  Let's decode the error flag into something more readable.

```@example tutorial_00_NLEVP
msg = decode_error_flag(flag)
println(msg)
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

```@example tutorial_00_NLEVP
sol,n,flag=mslp(T,0,tol=1E-10,output=true);
nothing #hide
```

Boom! The iteration stops after 6 iterations and indicates that an
eigenvalue has been found.

We may also inspect the corresponding left and right eigenvector. They are
fields in the solution object:

```@example tutorial_00_NLEVP
println(sol.v)
println(sol.v_adj)
```

There are more iterative solvers. Not all of them solve the complete problem.
Some only return the eigenvalue and the right eigenvector others even just the
eigenvalue. (If the authors of this software get more funding this might be
improved in the future...)

Let's briefly test them all...

### inverse iteration

This algorithm finds an eigenvalue and a corresponding right eigenvector

```@example tutorial_00_NLEVP
sol,n,flag=inveriter(T, 0, tol=1E-10, output=true);
nothing #hide
```

### trace iteration
This algorithm only finds an eigenvalue

```@example tutorial_00_NLEVP
sol,n,flag=traceiter(T, 0, tol=1E-10, output=true);
nothing #hide
```

### Lancaster's generalized Rayleigh quotient iteration
This algorithm only finds an eigenvalue

```@example tutorial_00_NLEVP
sol,n,flag=lancaster(T, 0, tol=1E-10, output=true);
nothing #hide
```

### two-sided Rayleigh-functional iteration
This algorithm finds a complete eigentriple.

```@example tutorial_00_NLEVP
sol,n,flag=rf2s(T, 0, tol=1E-10, output=true);
nothing #hide
```

Oops, this last test produced a NaN value. The problem often occurs when
the eigenvalue and one point in the iteration is too precise. (See how in the
last-but-one iteration we've been already damn close to 1/3). Choosing a
slightly different initial guess often leverages the problem:

```@example tutorial_00_NLEVP
sol,n,flag=rf2s(T, 0, tol=1E-10, output=true);
nothing #hide
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

```@example tutorial_00_NLEVP
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

```@example tutorial_00_NLEVP
Λ,V=beyn(T,Γ,l=6, output=true);
nothing #hide
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

```@example tutorial_00_NLEVP
for (idx,λ) in enumerate(Λ)
   v=V[:,idx]
   v/=sqrt(v'*v)
   res=T(λ)*v
   println("λ=$λ residual norm: $(abs(sqrt(res'*res)))")
end
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

```@example tutorial_00_NLEVP
number=count_poles_and_zeros(T,Γ)
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

