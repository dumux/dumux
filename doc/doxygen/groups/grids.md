# Doxygen subgroup for grids in the input output group

@defgroup Grids Grids
@brief Grids available in DuMu<sup>X</sup>
@ingroup InputOutput

DuMu<sup>X</sup> interacts with grids via the DUNE grid interface \cite BASTIAN2008. It supports grid implementations that implement this interface. Below you find a non-exhaustive list of several such grids from the DUNE community. (Typically, using grids requires \ref external-dependencies ["installing the respective grid module"].)

For the unified construction of grids, DuMu<sup>X</sup> provides adapters for several grids. These are specializations of the Dumux::GridManager template for the specific Grid type. You find these classes listed below. Each grid manager class's documentation lists the necessary and optional parameter keys to set in the input parameter file.

The following table provides an overview of grids with an existing Dumux::GridManager specialization and their capabilities.

| Grid | Type | Create grid directly from input file | Create grid from grid file | Grading | Periodic | Overlap | Partitioning | Refinement |
|------|------|--------------------------------------|----------------------------|----------|----------|----------|--------------|------------|
| \ref yasp_grid_manager "YASP" | Dune::YaspGrid<dim, Coordinates> | yes | dgf/vti | no | yes | yes | yes | yes |
|      | Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ctype, dim>> | yes | no | yes | yes | yes | yes | yes |
| \ref oned_grid_manager "OneDGrid" | Dune::OneDGrid | yes | dgf | no | no | no | no | yes |
| \ref ug_grid_manager "UGGrid" | Dune::UGGrid<dim> | yes | dgf/msh/vtp/vtu | no | no | no | no | yes |
| \ref alu_grid_manager "ALUGrid" | Dune::ALUGrid<dim, dimworld, elType, refinementType> | yes | dgf/msh/vtp/vtu | no | no | no | no | yes |
| \ref foam_grid_manager "FoamGrid" | Dune::FoamGrid<dim, dimworld> | yes | dgf/msh/vtp/vtu | no | no | no | no | yes |
| \ref porenetwork_grid_manager "Pore Network" ¹ | Dune::FoamGrid<dim, dimworld> | yes | dgf | no | no | no | no | no |
| \ref sp_grid_manager "SPGrid" | Dune::SPGrid<ct, dim, Ref, Comm> | yes | dgf | no | yes | yes | no | yes |
| \ref mmesh_grid_manager "MMesh" | MMesh | yes | dgf/msh/vtp/vtu | no | no | no | no | yes |

Additionally, \ref sub_grid_manager "dune-subgrid" can be used with any kind of host grid, allowing a subset of its elements to be used for the simulation.


¹ The pore network grid is based on FoamGrid and can be generated from either an input file or a DGF file. Parameters for throats and pore bodies may be defined globally in the input file or individually within the DGF file. It is also possible to generate a structured grid with random connections, where the diameters follow a specified distribution.

@defgroup PoreNetworkModelGrids Pore network model grids
@ingroup Grids
