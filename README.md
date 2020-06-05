# WavesAndEigenvalues.jl
*Julia package for handling various wave-equations and (non-linear) eigenvalue problems.*

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JulHoltzDevelopers/WavesAndEigenvalues.jl/dev)

## Disclaimer
The current version is not final and there are still some features to add, clean-up to be done, and open documentation tasks.

However, it's a start. So have fun!

## What is WavesAndEigenvalues.jl?
The package has evolved from academic research in thermoacoustic-stability analysis. But it is designed fairly general.
The package currently contains three modules, each of which is targeting at one of the main design goals:

1. Provide an elaborate interface to define, solve, and perturb nonlinear eigenvalues (**NLEVP**)
2. Provide a lightweight interface to read unstructured tetrahedral meshes in gmsh's version 4 file format. (**Meshutils**)
3. Provide a convenient interface for solving the (thermo-acoustic) Helmholtz equation. (**Helmholtz**)

## NLEVP
Assume you are a *literally* a rocket scientist and you want to solve an eigenvalue problem like
```
(K+omega*C + omega^2*M+ n*exp(-im*omega*tau)*F]p = 0
```
Where `K`,`M`, and `F` are some matrices, `im` the imaginary unit, `n` and `tau` some parameters, and `omega` and `p` unknown eigenpairs.
Then the **NLEVP** module lets you solve this equation and much more. Actually, any non-linear eigenvalue problem can be solved. Doing science in some other field then rockets or being pure mathematician is, thus, no problem at all. Just specify your favourite nonlinear eigenvalue problem and have fun. And *fun* here includes adjoint-based perturbation of your solution up to 30th order!  

## Meshutils
OK, this was theoretical. But you are a real scientist, so you know that the matrices `K`,`C`,`M` and `F` are obtained by discretizing some equation using a specific mesh. You do not want to be too limited to simple geometries so unstructured tetrahedral meshing is the method of choice.  **Meshutils** is the module that gives you the tools to read and process such meshes.

## Helmholtz
You are a practitioner? Fine, then eigenvalue theory and mesh handling should not bother your every day work too much.
You just want to read in a mesh, specify some properties (boundary conditions, speed of sound field...), and get a report on
the stability of your configuration? Then **Helmholtz** is your module! It even tells you how to optimize your design.
