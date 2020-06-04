# WavesAndEigenvalues.jl

*Documentation for WavesAndEigenvalues.*

A package for building Finite Element Models of wave-based problems and solving the associated eigenvalue problems.

!!! note

        This Package is under development, we will add more features soon!

## Package Features

- Load mesh files in .msh format, with surface and volumes labeled (see e.g. [gmsh](https://gmsh.info/)).
- Specify the equations to be solved on each surface and volumes.
- Build a paramteric depedent FEM `sparse_matrix`.
- Calculate eigenvalues of nonlinear eigenvalue problems (NLEVP).

The **Manual** documents the package modules and functionalities.

Examples on usage in the form of Jupyter Notebooks will be made available soon.

## Documentation Outline

```@contents
Pages = [
        "Mesh.md",
        "Helmholtz.md",
        "NLEVP.md",
]
Depth = 2
```

## Examples

```@contents
Pages =[
        "load_mesh.md",
]
Depth = 1
```
