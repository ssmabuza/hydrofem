
# Introduction

The "Hydrofem" project is a flexible serial finite element solver based on scalar nodal continuous finite elements. The code is designed to be a learning
tool and uses third party libraries `Eigen` and `CLI`. 

# Third Party Libraries (TPLs)

`Eigen` is already included in this repo. To get `CLI`, run either one of the following command inside the `tpls` directory:

 - `git clone git@github.com:CLIUtils/CLI11.git`
 - `git clone https://github.com/CLIUtils/CLI11.git`

# Build Instructions

After cloning this repository, create a build directory, say `build`. In the build directory run the following:

 - `cmake ..`
 - `make -j2`

# Running Simulations

Once built, the executable `hydrofem` maybe be used to get various simulations started depending on the options provided either on the command line or in a configuration script. Here are some examples:

 - In general a script would be run this way: `/path/to/hydrofem  --config=/path/to/config/file`
 - Specific example if in the build directory: `./hydrofem --config=config.poisson`
