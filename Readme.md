
# Introduction

The "Hydrofem" project is a flexible partly serial and partly parallel finite element solver based on scalar nodal continuous finite elements. The code is designed to be a learning
tool and uses third party libraries: `Boost`, `Eigen`. The project will swap out the `Eigen` library solvers and use the parallel `PETSc` and `mfem` libraries. 

# Third Party Libraries (TPLs)

Install `Boost` via the linux package managers or compile it from source. In Ubuntu/Debian, write `sudo apt-get install libboost-all-dev`. `Eigen` is already included in this repo. 
`mfem` and `PETSc` are installed together using the [spack](https://spack.io/) framework.

# Build Instructions

After cloning this repository, use cmake presets to generate the makefiles and compile the project by running:

 - `cmake --preset debug` for the debug mode. Replace debug with release for the release mode.
 - `make --build --preset debug` to compile the project.

# Running Unit Tests

The unit test executable is found in the build directory.

# Running Simulations

Once built, the executable `hydrofem` maybe be used to get various simulations started depending on the options provided either on the command line or in a configuration script. Here are some examples:

 - In general a script would be run this way: `/path/to/hydrofem  --config=/path/to/config/file`
 - Specific example if in the build directory: `./hydrofem --config=config.poisson`
