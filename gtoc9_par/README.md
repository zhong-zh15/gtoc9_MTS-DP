MultiTree Search is Open Source algorithm for Granular Dynamics Simulations for Astrodynamics and Planetary Science. 

DEMBody stands for Discrete Element Model Body, which incorporates classical granular dynamics simulations and N-body self-gravity calculation. The software focus on investigation of geological features, surface evolution and in-situ exploration on small celestial bodies, e.g., size sorting/segregation on asteroids, mass creeping/wasting driven by geological processes and locomotion dynamics of the small body lander/rover on small body granular regolith. And this code can be easily customized to accurately and efficiently solve other problems in astrophysics, e.g., simulating planetesimals, moons, ring or dust particles.

The software implements the Soft-Sphere Discrete Element Model (SSDEM), coupled with the N-body gravity integrator. Parallized with OpenMP (for shared memory systems) and MPI (for distributed memory systems), DEMBody can execute on supercomputers, clusters or multi-core PCs running a Linux-based (or Windows, but not recommended) operating system. We developed several branches for DEMBody with various special features, e.g., DEMBody-bond with bonded spheres, DEMBody-stretchPBC with moving periodic boundary conditions and DEMBody-MixBiMesh with two lattice structures with different cell sizes for particle systems with a bidisperse distribution.

If you use this code or parts of this code for results presented in a scientific publication, we would greatly appreciate a citation eithor to this code or to our published papers listed below.

## Papers based on DEMBody

* [Reconstructing the formation history of top-shaped asteroids from the surface boulder distribution](https://doi.org/10.1038/s41550-020-01226-7) - *Nature Astronomy*, 2021
* [Numerical simulations of the controlled motion of a hopping asteroid lander on the regolith surface](https://doi.org/10.1093/mnras/stz633) - *MNRAS*, 2019
* [Collision-based understanding of the force law in granular impact dynamics](https://doi.org/10.1103/PhysRevE.98.012901) - *PRE*, 2018
* [Asteroid surface impact sampling: dependence of the cavity morphology and collected mass on projectile shape](https://doi.org/10.1038/s41598-017-10681-8) - *Scientific Reports*, 2017

## Getting Started

These instructions will get you a copy of the project for development and testing purposes.

### Structure
  `Data`: storage for data file including point data, wall data and other miscellaneous data (bond data, biDisperse data).

  `Input`: input files including systemControl.dembody, input_points.txt, bondedWallPoint.vtk, trimeshWall.mesh, bondedTriMeshWall.mesh, largeParticles.bidisperse, gravTriMesh.force.

  `src`: code of DEMBody.

  `build`:  build DEMBody by Cmake.

  `example`: a series of examples showing some features in DEMBody.

  `doc`: instruction on basic models and tutorial for the code.

### Quick installation 

You can run DEMBody by only steps:

```
git clone --recursive https://github.com/Bin-Cheng-THU/DEMBody
cd DEMBody/build && rm -rf *
cmake ../src
make
cp ../src/DEMBody.sh . && sbatch DEMBody.sh
```

### Documentation 

The online documentation with many examples and tutorials can be found at [readthedocs.io](dembody.readthedocs.io).

### Pre-/Post-process

We developed [DEMTool](https://github.com/Bin-Cheng-THU/DEMTool.git) for pre- and post-process for DEMBody, e.g., point/mesh file generation and data rending based on [POV-Ray](http://www.povray.org/). It is a standalone package whose output files are easily imported into most commercial and Open Source DEM software. Please read the instruction in [DEMTool](https://github.com/Bin-Cheng-THU/DEMTool.git) for details.

## Authors

* **Bin Cheng** - *Initial work* - [Bin-Cheng-THU](https://github.com/Bin-Cheng-THU)

See also the list of [contributors](AUTHORS.md) who participated in this project.

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE - see the [LICENSE](LICENSE) file for details

<!---
## Acknowledgments

* Prof. Baoyin and colleagues in LAD
* My girlfriend Fanbing Zeng
* My parents, brother and whole family
-->