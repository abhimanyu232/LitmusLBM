```
██╗     ██╗████████╗███╗   ███╗██╗   ██╗███████╗      ██╗     ██████╗ ███╗   ███╗
██║     ██║╚══██╔══╝████╗ ████║██║   ██║██╔════╝██╗██╗██║     ██╔══██╗████╗ ████║
██║     ██║   ██║   ██╔████╔██║██║   ██║███████╗╚═╝╚═╝██║     ██████╔╝██╔████╔██║
██║     ██║   ██║   ██║╚██╔╝██║██║   ██║╚════██║██╗██╗██║     ██╔══██╗██║╚██╔╝██║
███████╗██║   ██║   ██║ ╚═╝ ██║╚██████╔╝███████║╚═╝╚═╝███████╗██████╔╝██║ ╚═╝ ██║
╚══════╝╚═╝   ╚═╝   ╚═╝     ╚═╝ ╚═════╝ ╚══════╝      ╚══════╝╚═════╝ ╚═╝     ╚═╝
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░                                                    
```
## A High Performance C++ framework for Computational Fluid Dynamics using the Lattice Boltzmann Method

currently under active development.  

**Dependencies and Requirements**
------
* Compiler : Tested with `gcc 15.2.0` and above, with `--std=c++23`

**Current Model Features**
---
- Standard Incompressible LBGK
  - Lattices
    - 2D: D2Q5, D2Q9
    - 3D: D3Q19, D3Q27 
  - Validated Setups
      - Taylor Green Vortex 2D
      - Double Periodic Shear Layer 2D

**Build instructions** 
---
1. Configure CMake: `cmake -B build -DCMAKE_BUILD_TYPE={Release|RelWithDebInfo|Debug}`
2. Build and Compile: ` cmake --build build `
    * Default CMAKE_BUILD_TYPE = RelWithDebInfo


**Upcoming Features and Upgrades**
---
- 3D Validation
- Boundary Conditions
- Shared Memory Parallelisation (OpenMP)
- Compressible Flow Simulations with Double Distribution Function Models
- Cuda Implementation


**Author**
---
Dr. Abhimanyu Bhadauria, Ph.D
