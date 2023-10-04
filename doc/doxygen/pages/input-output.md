# Input and output formats

This section explains summarizes
the grid formats that can be used by DuMux
and introduces grid management via Dumux::GridManager.
We briefly discuss grid generation. Finally, we discuss the options for
outputting results for checkpointing or visualization.

## Supported grid formats

DuMux can read grids from file using the Dune Grid Format (DGF),
the Gmsh mesh format (MSH, version 2),
the Eclipse grid format (GRDECL),
and various Visualization ToolKit (VTK/VTU/VTP) formats.

### The Dune Grid Format (DGF)

Most of our DuMux tests use the Dune Grid Format (DGF) to read in grids. A detailed description
of the DGF format and some examples can be found in the DUNE doxygen documentation
([Modules $\rightarrow$ I/O $\rightarrow$ Dune Grid Format (DGF)](https://www.dune-project.org/doxygen/2.9.0/group__DuneGridFormatParser.html#details)). To generate larger or more
complex DGF files, we recommend to write your own scripts, e.g, in C++, Matlab or Python.

The DGF format can also be used to read in spatial parameters defined on the grid. These parameters can
be defined on nodes as well as on the elements. An example for predefined parameters on a grid
can be found in [`dumux/test/porousmediumflow/co2`](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/tree/master/test/porousmediumflow/co2).

### The Gmsh MSH format

[Gmsh](http://geuz.org/gmsh/) is a powerful open-source grid generator for unstructured (and structured)
finite-element meshes (@cite GEUZAINE2009).
DuMux supports the default Gmsh mesh format (MSH), version 2, through Dune.
For the format specifics and how to create grids with Gmsh, e.g., using
the provided GUI, we refer to the [Gmsh documentation](http://geuz.org/gmsh/doc/texinfo/gmsh.html).
Note that the current default in Gmsh produces files in the format version 4. Use the option

```sh
gmsh -format msh2
```

or add `Mesh.MshFileVersion = 2.2;` in a `.geo` file.
The MSH format can contain element and boundary markers defined on the grid.
Thus, boundaries can be easily marked as, e.g., inflow boundaries
using Gmsh. Further, the format supports higher order elements.
They can be used to create boundary parametrization supported by, e.g., the grid
managers `Dune::UGGrid` and `Dune::FoamGrid`.

### Eclipse Grid Format

The Eclipse Grid Format (GRDECL) is commonly used for corner-point grids.
Such grids consist of hexahedra, which are described by eight points on so-called pillars.
A special feature of corner-point geometries is that points on pillars can degenerate,
meaning that two neighboring points on a pillar can coincide.
Furthermore, faces are, in general, bi-linear and cells can be non-convex.
This allows for the accurate description of faults, layers, or wells, occurring in geological environments.

Furthermore, petrophysical properties can be defined (for each cell),
by using eclipse-specific keywords, e.g. `PORO`, `PERMX`, `PERMY`.

DuMux supports the Eclipse Grid Format by using the [`opm-grid` module](https://github.com/OPM/opm-grid)
of the [Open Porous Media initiave (OPM)](https://opm-project.org).
See @ref external-libraries for how to install `opm-grid` together with DuMux.
An example using a corner-point grid can be found in `dumux/test/porousmediumflow/2p/cornerpoint`.

### VTK File Format

VTK format uses ASCII or XML format. It is mostly used by DuMux for output purposes and can be visualized
by programs such as Paraview, VisIt or Tecplot. Using VTK files to input grid and parameter data is also possible.
An example can be found in `dumux/test/io/gridmanager`.

### Other Grid Formats

Grid formats other than DGF, MSH, GRDECL, or VTK will have to be converted to the DGF, MSH, GRDECL, or VTK format before they can be used in DuMux.
A powerful converter tool supporting many common formats is the Python tool [`meshio`](https://pypi.org/project/meshio/).
If conversion is not an option, another possibility would be to write your own manager class based on Dumux::GridManager.

## The DuMux Grid manager class

The Dumux::GridManager class helps to construct the grid from information in the input file and also handles data provided in grid files.
Currently, supported Dune grid interface implementations are `Dune::YaspGrid`, `Dune::OneDGrid`, `Dune::UGGrid`, `Dune::ALUGrid`, `Dune::FoamGrid`, `Dune::SubGrid`, `OPM::Grid` (cornerpoint grids), and `Dune::SPGrid`
(also see [overview of grid managers in Dune](https://www.dune-project.org/groups/grid/)).
Grids can be constructed from a DGF, VTK or MSH file by simply providing the filename to the grid in the `[Grid]` group.
Note that group name `[Grid]` is the default group name and can be customized in your problem changing the string property Dumux::Properties::GridParameterGroup. For setups with more than one grid, a different group name can be set for each grid,
which allows to set different options per grid.

An exemplary input file configuring the grid manager to read a grid from the file `mydgfgrid.dgf` looks like this

```ini
[Grid]
File = mydgfgrid.dgf
```

If you are using an unstructured grid interface like `UGGrid` or `FOAMGrid`, constructing a grid from a VTK or MSH is just changing a line:

```ini
[Grid]
File = mygmshgrid.msh
```

DuMux will tell you in case your selected grid manager does not support reading such files.

You want to initially refine your grid? It's just adding a line:

```ini
[Grid]
File = mydgfgrid.dgf
Refinement = 4
```

When reading a MSH or VTK file, further parameters are recognized. `Verbose` enables verbose output on grid construction when set to $1$.
`BoundarySegments` enables reading parameterized boundaries. `PhysicalEntities` enables reading boundary and element flags.

### Parameters specific to the grid implementation

The Dumux::GridManager supports also a selection of parameters that are specific to the chosen grid manager.
To give an example, we take a look at the unstructured grid `Dune::UGGrid`.
`UGGrid` supports red-green refinement per default. One can turn off the green closure by setting the grid's closure type

```ini
[Grid]
File = mydgfgrid.dgf
ClosureType = None # or Green
```

For all available parameters, see Dumux::GridManager and its specializations for different grid types.

### Structured grids
If you want to construct a structured grid without using a specific grid file, insert the following into the input file:

```ini
[Grid]
LowerLeft = 0 0 0
UpperRight = 1 1 1
Cells = 10 10 20
```

where `LowerLeft` is a vector to the lower left corner of the grid and `UpperRight` a vector to the upper right corner.
`Cells` is a vector with the number of cells in each coordinate direction. Note,  that for a grid in a two-dimensional world, the
vectors only have two entries.

Depending on the grid manager, further parameters are recognized.
`UGGrid`, for example, supports simplex elements as well as hexahedral elements
(called ``cube'' in Dune). When creating a structured grid, we can select the cell type as follows

```ini
[Grid]
LowerLeft = 0 0 0
UpperRight = 1 1 1
Cells = 10 10 20
CellType = Cube # or Simplex
```

For all available parameters see Dumux::GridManager<Dune::YaspGrid<dim,Coordinates>>
and Dumux::GridManager<Dune::YaspGrid<dim,Dune::TensorProductCoordinates<ctype,dim>>> class documentations.

### Other grid manager implementations

* Dumux::CakeGridManager: Provides a method to create a piece of cake grid
* Dumux::CpGridManager: Reads the GRDECL file and generates a corner-point grid

## Other input and output formats

The following formats are supported for visualization, and general
data input and output.

### VTK file format (output)

Dumux allows to write out simulation results via the `VtkOutputModule`.
For every print-out step, a single VTU file is created. For parallel simulations one file
per print-out step is generated for each processor.
The PVD file groups the single VTU files and contains additionally the time step information.
The VTK file format is supported by common visualisation programs like ParaView, VisIt, and Tecplot.

If provided by the model implementation, the `initOutputModule` function of the model's `IOFields`,
adds a default set of variables to the VTK output module instance. It is also possible to add variables,
using the member function Dumux::VtkOutputModuleBase::addField of the Dumux::VtkOutputModule.
For example, to add a variable called `temperatureExact` add

```cpp
vtkWriter.addField(problem->getExactTemperature(), "temperatureExact");
```

The first input argument of this method is the value of the additional variable,
provided by a method of the corresponding problem.
If it does not already exists, the user has to provide this method.

```cpp
const std::vector<Scalar>& getExactTemperature()
{ return temperatureExact_; }
```

It is important that the life-time of the added field exceeds the life-time of the writer. That means you cannot pass temporaries
to the Dumux::VtkOutputModuleBase::addField function. The vector has to be stored somewhere, e.g. in the program main function.

The second input argument is the name of the additional variable (as it should be written in the VTK files).
The example above is taken from `test/porousmediumflow/1pnc/implicit/1p2c/nonisothermal/convection/main.cc`

### VTK file format (input)
There is support for reading data and grids from VTK files, see Dumux::VTKReader.

## Gnuplot interface
DuMux provides a small interface to [Gnuplot](http://www.gnuplot.info/),
which can be used to plot results and generate
image files (e.g. `.png`). To use the gnuplot, gnuplot has to be installed.
The following is a brief introduction. For the class documentation, see Dumux::GnuplotInterface.

A Gnuplot interface is available to plot or visualize results during a simulation run.
To use the gnuplot interface you have to make some modifications in your file, e.g., your main file.

First, you have to include the corresponding header file for the gnuplot interface.

```cpp
#include <dumux/io/gnuplotinterface.hh
```

Second, you have to create an instance of the class
Dumux::GnuplotInterface (e.g. called `gnuplot`).

```cpp
Dumux::GnuplotInterface<double> gnuplot;
```

As an example, to plot the mole fraction of nitrogen (`y`) over time (`x`),
extract the variables after each time step in the time loop.
The actual plotting is done using the method of the Gnuplot interface:

```cpp
gnuplot.resetPlot();                         // reset the plot
gnuplot.setXRange(0.0, 72000.0);             // specify xmin and xmax
gnuplot.setYRange(0.0, 1.0);                 // specify ymin and ymax
gnuplot.setXlabel("time [s]");               // set xlabel
gnuplot.setYlabel("mole fraction mol/mol");  // set ylabel
// set x-values, y-values, the name of the data file and the Gnuplot options
gnuplot.addDataSetToPlot(x, y, "N2.dat", options);
gnuplot.plot("mole_fraction_N2");            // set the name of the output file
```

It is also possible to add several data sets to one plot by calling Dumux::GnuplotInterface::addDataSetToPlot more than once.
For more information have a look into a test including the gnuplot interface header file, the class documentation
of Dumux::GnuplotInterface, or the header file itself `dumux/io/gnuplotinterface.hh`.

## Container I/O
DuMux supports writing to file from and reading
into some STL containers like `std::vector<double>` or `std::vector<Dune::FieldVector>`.
If you want to read and write simple vectors, have a look at the header `dumux/io/container.hh`.

## Matrix and Vector I/O

`dune-istl` (the Dune Iterative Solver Template Library) supports writing and reading vectors
and matrices to/from different format. For example you can write a matrix in a sparse matrix format that
can be read by Matlab (see the header `dune/istl/io.hh`).

## Restarting simulations (check-pointing)

DuMux has some experimental support for check-pointing (restarting paused/stopped/crashed simulations).
You can restart a DuMux simulation from any time point where a VTK file was written out.
This is currently only supported for sequential, non-adaptive simulations. For adaptive simulation
the full hierarchical grid has to be stored. This is usually done with the grid's `Dune::BackupRestoreFacility`.
There is currently no special support by DuMux for that, but it is possible to implement
a restart using `Dune::BackupRestoreFacility` with plain Dune.

For VTK files the output can be read with the free function `Dumux::loadSolution`. Grids can be read with
the `Dumux::VTKReader` or you can simply recreate the grid as you did in the first simulation run.

Writing double-precision floating point numbers to VTK files is available since Dune release 2.7.
If you are using that version, it is now possible to specify output precision in the input file using
`Dumux::Vtk::Precision` followed by either `Float32`, `Float64`, `UInt32`, `UInt8` or `Int32`
`Float32` is set as the default. We especially advice the use of `Float64` when working with restart files.

The restart capabilities will hopefully be improved in future versions of DuMux $3$.
We are looking forward to any contributions (especially HDF5 / XDMF support, improvement of VTK support).
