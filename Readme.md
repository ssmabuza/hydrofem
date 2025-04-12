
# Introduction

The "Hydrofem" project is a flexible partly serial and partly parallel finite element solver based on scalar nodal continuous finite elements. The code is designed to be a learning
tool and uses third party libraries: `Boost`, `Eigen`, `CLI`, `PETSc` and `deal.ii`. 

# Third Party Libraries (TPLs)

Install `Boost` via the linux package managers or compile it from source. In Ubuntu/Debian, write `sudo apt-get install libboost-all-dev`. `Eigen` is already included in this repo. To get `CLI`, run either one of the following command inside the `tpls` directory:

 - `git clone git@github.com:CLIUtils/CLI11.git`
 - `git clone https://github.com/CLIUtils/CLI11.git`

`deal.ii` and `PETSc` are installed together using the [spack](https://spack.io/) framework.

# Build Instructions

After cloning this repository, create a build directory, say `build`. In the build directory run the following:

 - `cmake ..`
 - `make -j2`

# Running Simulations

Once built, the executable `hydrofem` maybe be used to get various simulations started depending on the options provided either on the command line or in a configuration script. Here are some examples:

 - In general a script would be run this way: `/path/to/hydrofem  --config=/path/to/config/file`
 - Specific example if in the build directory: `./hydrofem --config=config.poisson`
