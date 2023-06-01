# Directory Structure

DuMux has the following folder structure, which is similar to other DUNE modules.

* `bin`: binaries, e.g. used for the automatic testing, post-processing, installation
* `cmake`: the configuration options for building DuMux
* `doc`: files necessary for the Doxygen documentation and various logos
* `dumux`: the main folder, containing the source files. See below for a visualized structure.
* `examples`: well-documented examples of applying DuMux to typical physical problems. In the `README.md` files, the setup is explained, the used code is presented as well as documented and images resulting from the simulation are included. The `README.md` files are located directly in the subfolders of `examples` and can be displayed by web browsers.
* `test`: tests for each numerical model and some functionality.The structure is equivalent to the `dumux` folder, the `references` folder contains solutions for the automatic testing. Each test program consist of a main file `main.cc`, the problem definition `*problem.hh` (specifying initial and boundary conditions), and an input file `params.input`. If necessary, spatially dependent parameters are defined in `*spatialparameters.hh`. For more detailed descriptions of the tests, please have a look at the [examples](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/tree/master/examples#open_file_folder-example-1-diffusion-equation).

@mermaid{"./images/dumux_structure.mmd}