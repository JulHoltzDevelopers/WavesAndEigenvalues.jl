```@meta
EditURL = "<unknown>/tutorial_07_Bloch_periodicity.jl"
```

# Tutorial 07 Bloch periodicity

The Helmholtz solver was designed with the application of thermoacoustic
stability assesment in mind. A very common case is that an annular combustion
chamber needs assessment. These chamber types feature a discrete rotational
symmetry, i.e., they are built from one cionstituent entity -- the unit cell -- that
implicitly defines the annular geometry by copying it around the circumference
multiple times. This high regularity of the geometry can be used to efficiently
model the problem, e.g. only storing the unit cell instead of the
entire geometry in a mesh-file in order to save memory and even performing the
entire analysis just and thereby significantly accelerating the computations.[1]
This tutorial teaches you how to use these features.

!!! note
    To do this tutorial yourself you will need the `"NTNU_12.msh`" file.
    Download it [here](NTNU_12.msh).

## Model Set-up

As a geometry we choose the laboratory-scale annular combustor used at NTNU
Norway [2]. The NTNU combustor features unit cells that are reflection-symmetric.
The meshfile `NTNU_12.msh`, thus, stores only one half of the unit cell. Let's
load it...

```@example tutorial_07_Bloch_periodicity
using WavesAndEigenvalues.Helmholtz
mesh=Mesh("NTNU_12.msh",scale=1.0);
nothing #hide
```

Note that there are two special domains named `"Bloch"` and `"Symmetry"`. These
are surfaces associated with the special symmetries of the combustor geometry.
`"Bloch"` is the surface where one unit cell ends and the next one should begin,
while `"Symmetry"` is the symmetry plane of the unit-cell. As the mesh is stored
as a half-cell both surfaces are boundaries. The code will automatically use
them to create the full mesh or just a unit-cell. First ,let's create a full mesh
We therefore specify the expansion scheme, this is a list of tuples containing
all domain names we like to expand together with a qualifyer that is either
`:full`,`:unit`, or `:half`, in order to tell the the code whether the copied
versions should be merged into one addressable domain of the same name spanning
the entire annulus (`:full`), pairs of domains from two corresponding half cells
should be merged into individually adressable entities (`:unit`), or the copied
hald cells remain individually addressable (:half). In the current case we only
need individual names for each flames and can sefely merge all other entities.
The scheme therefore reads:

```@example tutorial_07_Bloch_periodicity
doms=[("Interior",:full),("Inlet",:full), ("Outlet_high",:full),
     ("Outlet_low",:full), ("Flame",:unit),];
nothing #hide
```

and we can generate the mesh by

```@example tutorial_07_Bloch_periodicity
full_mesh=extend_mesh(mesh,doms);
nothing #hide
```

The code automatically recognizes the degree of symmetry -- 12 in the current
case-- and expands the mesh. Note how the full mesh now has 12 connsecutively
numbered flames, while all other entities still have a single name even though
they were copied expanded.

## Discretizing the full mesh

The generated mesh behaves just like any other mesh. We can use it to discretize
the model.

```@example tutorial_07_Bloch_periodicity
speedofsound(x,y,z)=z<0.415 ? 347.0 : 850.0; #speed of sound field in m/s
C=generate_field(full_mesh,speedofsound);
dscrp=Dict(); #model description
dscrp["Interior"]=(:interior,());
dscrp["Outlet_high"]=(:admittance, (:Y_in,0));# src
dscrp["Outlet_low"]=(:admittance, (:Y_out,0));# src
L=discretize(full_mesh, dscrp, C,order=:lin)
```

We can use all our standard tools to solve the model. For instance, Indlekofer
et al. report on a plenum-dominant 1124Hz-mode that is first order with respect to the
azimuthal direction and also first order with respect to the axial direction.
Let's see whether we can find it by putting 1000 Hz as an initial guess.

```@example tutorial_07_Bloch_periodicity
sol,nn,flag=mslp(L,1000,tol=1E-9,scale=2pi,output=true);
nothing #hide
```

There it is! But to be honest the computation is kinda slow. Obviously, this is
because the full mesh is quite large. Let's see how the unit cell computation
would do....

## Unit cell computation

Ok, first, we create the unit cell mesh. We therefore set the keyowrd parameter
`unit` to `true`in the call to `extend_mesh`.

```@example tutorial_07_Bloch_periodicity
unit_mesh=extend_mesh(mesh,doms,unit=true)
```

Voila, there we have it! Discretization is nearly as simple as with a normal mesh.
However, to invoke special Bloch-periodic boundary conditions the optional
parameter `b` must be set to define a symbol that is used as the Bloch wavenumber.
This option will also triger. The rest is as before.

```@example tutorial_07_Bloch_periodicity
#
c=generate_field(unit_mesh,speedofsound);
l=discretize(unit_mesh, dscrp, c, b=:b)
```

Note how the signature of the discretized operator is now much longer. These
extra terms facilitate the Bloch periodicity. What is important is the parameter
`:b` that is used to set the Bloch wavenumber. For relatively low frequencies,
the Bloch wave number corresponds to the azimuthal mode order. Effectively, it
acts like a filter to the eigenmodes. Setting it to `1` will give us mainly first
order modes.   (In general it gives us modes with azimuthal order `k*12+1` where
`k` is a natural number including 0, but 13th order modes are very high in frequency
and therefore don't bother us). In comparison to the full model, this is an extra
advantage of Bloch wave theory, as we cannot filter for certain mode shapes using
the full mesh. This feature also allows us to reduce the estimate for eigenvalues
inside the integration contour when using Beyn's algorithm, as there can only be
eigenmodes corresponding to the current Bloch wavenumber. Let's try finding
1124-Hz mode again. It's of first azimuthal order so we should put the Bloch
wavenumber to 1.

```@example tutorial_07_Bloch_periodicity
l.params[:b] = 1;
sol,nn,flag = mslp(l,1000,tol=1E-9, scale=2pi, output=true);
nothing #hide
```

The computation is much faster than with the full mesh, but the eigenfrequency
is exactly the same!

## Writing output to paraview
Not only the eigemnvalues much but also the mode shapes. Bloch wave theory
clearly dictates how to extzend the mode shape from the unit cell to recover
the full-annulus solution.

The vtk_write function does this for you. You, therefore, have to provide it
with the current Bloch wavenumber. If you provide it with the unit_cell mesh,
it will just write a single sector to the file.

```@example tutorial_07_Bloch_periodicity
data=Dict();
v=bloch_expand(full_mesh,sol,:b);
data["abs"]=abs.(v)./maximum(abs.(v));
data["pahase"]=angle.(v)./pi;
vtk_write("Bloch_tutorial",full_mesh,data);
nothing #hide
```

## Summary

Bloch-wave-based analysis, significantly reduces the size of your problem
without sacrificing any precission. It thereby saves you both memory and
computational time. All you need to provide is

## References

[1] Efficient Computation of Thermoacoustic Modes in Industrial Annular
Combustion Chambers Based on Bloch-Wave Theory, G.A. Mensah, G. Campa,  J.P. Moeck
2016, J. Eng. Gas Turb. and Power,[doi:10.1115/1.4032335](https://doi.org/10.1115/1.4032335)

[2] The effect of dynamic operating conditions on the thermoacoustic response of
hydrogen rich flames in an annular combustor,  T. Indlekofer, A. Faure-Beaulieu,
N. Noiray and J. Dawson, 2021, Comb. and Flame, [doi:10.1016/j.combustflame.2020.10.013](https://doi.org/10.1016/j.combustflame.2020.10.013)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
