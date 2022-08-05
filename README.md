MultiTree Search is Open Source Algorithm for Multi-Spacecraft Successive Rendevous Mission. It was created in the 11th edition of the China Trajectory Optimization Competition (CTOC11) and then improved and tested for the 9th edition of the Global Trajectory Optimization Competition (GTOC9).

MultiTree Search is a natural extension of the original tree search framework, aiming to store with minimum space complexity, where the traditional tree search algorithm can be easily developed. This version combines the classic beam search with local search for better performance in GTOC9, which can be regarded as a variant multi-travel-salesman problem. We also apply dynamic programming to obtain the optimal rendevous time sequence (a continuous variable optimization problem). As far as we know, this is the first time we can give a theory-proof global optimal result in time optimization under a given rendevous sequence.

Parallized with OpenMP (for shared memory systems), MultiTree can execute on supercomputers, clusters, or multi-core PCs running a Linux-based (or Windows) operating system. Now, MultiTree can only be run on the C++ platform. 

If you use this code or parts of this code for results presented in a scientific publication, we would greatly appreciate a citation either to this code or to our published papers listed below.

## Papers based on Multi-Tree Search

* [Multi-Spacecraft Successive Rendevous: a Multi-Tree Search Embedded Dynamic Programming Algorithm
](https://doi.org/10.1038/s41550-020-01226-7) - *Journal of Guidance, Control, and Dynamics*, 2022 (Submitted)
* [Multi-Tree Search for Multi-Satellite Responsiveness Scheduling
Considering Orbital Maneuvering](https://doi.org/10.1093/mnras/stz633) - *IEEE TRANSACTIONS ON AEROSPACE AND ELECTRONIC SYSTEMS*, 2022

## Getting Started

These instructions will get you a copy of the project for development and testing purposes.

### Structure
  `input_data`: orbit elements of all space debris.
  
  `src`: code of MultiTree.

  `include`: h files needed.

  `output_result`: results computed by MultiTree

  `bin`: compiling file

  `build`:  build MultiTree by Cmake.

  `lib`: library needed (None in this version).

### Quick installation 

You can run MultiTree by only several steps:

* Linux 
```
git clone --recursive https://github.com/zhong-zh15/gtoc9_MTS-DP
cd gtoc9_par/build && rm -rf *
cmake ..
make
cd ../bin 
sbatch GTOC9_parallel.sh
```

* Windows (Visual Studio 2022)
```
git clone --recursive https://github.com/zhong-zh15/gtoc9_MTS-DP
```

## Authors

* **Zhong Zhang** - *Initial work* - [Zhong Zhang, Tsinghua, LAD](https://github.com/zhong-zh15)
* **Nan Zhang** - *Discussion and ideas verification*
* **Zherui Chen** - *Discussion and ideas verification*

<!---
See also the list of [contributors](AUTHORS.md) who participated in this project.
-->

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE - see the [LICENSE](LICENSE) file for details

<!---
## Acknowledgments

* Prof. Baoyin and colleagues in LAD
-->