# WavesAndEigenvalues.jl
Julia package for handling various wave-equations and (non-linear) eigenvalue problems.

## Installation

WavesAndEigenvalues can be installed from Julia's official package repository using the built-in package manager. Just type `]` to enter the package manager and then
```
pkg> add WavesAndEigenvalues
```
or
```julia-repl
julia> import Pkg; Pkg.add("WavesAndEigenvalues")
```
in the REPL and you are good to go.


## What is WavesAndEigenvalues.jl?
The package has evolved from academic research in thermoacoustic-stability analysis. But it is designed fairly general. It currently contains three modules, each of which is targeting at one of the main design goals:

1. Provide an elaborate interface to define, solve, and perturb nonlinear eigenvalues (**NLEVP**)
2. Provide a lightweight interface to read unstructured tetrahedral meshes using nastran (`*.bdf` and `*.nas`) or the latest gmsh format (`*.msh `). (**Meshutils**)
3. Provide a convenient interface for solving the (thermo-acoustic) Helmholtz equation. (**Helmholtz**)

## NLEVP
Assume you are a *literally* a rocket scientist and you want to solve an eigenvalue problem like
```math
(\mathbf K+\omega \mathbf C + \omega^2 \mathbf M+ n\exp(-i\omega\tau) \mathbf F]\mathbf p = 0
```
Where $\mathbf K$, $\mathbf C$, $\mathbf M$, and $\mathbf F$ are some matrices, $i$ the imaginary unit, $n$ and $\tau$ some parameters, and $\omega$ and $\mathbf p$ unknown eigenpairs.
Then the **NLEVP** module lets you solve this equation and much more. Actually, any non-linear eigenvalue problem can be solved. Doing science in some other field than rockets or even being a pure mathematician is, thus, no problem at all. Just specify your favourite nonlinear eigenvalue problem and have fun. And *fun* here includes adjoint-based perturbation of your solution up to arbitrary order!  

## Meshutils
OK, this was theoretical. But you are a real scientist, so you know that the matrices $\mathbf K$, $\mathbf C$, $\mathbf M$, and $\mathbf F$ are obtained by discretizing some equation using a specific mesh. You do not want to be too limited to simple geometries so unstructured tetrahedral meshing is the method of choice.  **Meshutils** is the module that gives you the tools to read and process such meshes.

## Helmholtz
You are a practitioner? Fine, then eigenvalue theory and mesh handling should not bother your every day work too much.
You just want to read in a mesh, specify some properties (boundary conditions, speed of sound field...), and get a report on
the stability of your configuration? **Helmholtz** is your module! It even tells you how to optimize your design.
