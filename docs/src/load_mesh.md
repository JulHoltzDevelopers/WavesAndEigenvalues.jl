# How to load a mesh?

First load tools from the Mesh utility module
    using WavesAndEigenvalues.Meshutils
Then load your mesh
    Mesh("filename.msh",scale=0.001)
The `scale` parameter is optional and multiplies the coordinates of the red mesh with some number. This functionality is basically here to convert units of length.
