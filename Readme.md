
# Introduction

The "Hydrofem" project is a flexible serial finite element solver based on scalar nodal continuous finite elements. The code is designed to be a learning
tool and uses third party libraries `Eigen` and `CLI`. 

# Third Party Libraries (TPLs)

`Eigen` is already included in this repo. To get `CLI`, run the following command inside the `tpls` directory:

 - `git clone git@github.com:CLIUtils/CLI11.git`

# Build Instructions

After cloning this repository, create a build directory, say `build`. In the build directory run the following:

 - `cmake ..`
 - `make -j2`


