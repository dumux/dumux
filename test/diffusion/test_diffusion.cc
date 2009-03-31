#include "config.h"
#include <iostream>
#include <iomanip>
//#ifdef HAVE_UG
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/sgrid.hh>
//#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include <dumux/material/matrixproperties.hh>
//
#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocity.hh"
//#include "dumux/diffusion/fe/fediffusion.hh"
//
//#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
//
#include "dumux/diffusion/problems/uniformproblem.hh"
#include "dumux/fractionalflow/variableclass.hh"

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=2;

        // create a grid object
        typedef double NumberType;
        typedef Dune::SGrid<dim,dim> GridType;

        Dune::FieldVector<GridType::ctype,dim> L(0);
        Dune::FieldVector<GridType::ctype,dim> R(300);
        Dune::FieldVector<int,dim> N(2);
        GridType grid(N,L,R);

        //Uniform mat;
        Dune::Uniform mat;

        Dune::HomogeneousSoil<GridType, NumberType> soil;
//        Dune::HeterogeneousSoil<GridType, NumberType> soil(grid, "permeab.dat", true);
//        printvector(std::cout, *(soil.permeability), "permeability", "row", 200, 1);
//        soil.permeability.vtkout("permeability", grid);

        Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, mat, mat);

        typedef Dune::VariableClass<GridType, NumberType> VC;

        double initsat = 1;

        VC variables(grid,initsat);

        Dune::UniformProblem<GridType, NumberType, VC> problem(variables, mat, mat, soil, materialLaw);

        Dune::Timer timer;
        timer.reset();
        //Dune::FEDiffusion<GridType, NumberType> diffusion(grid, problem);
        //Dune::FVDiffusion<GridType, NumberType, VC> diffusion(grid, problem);
        Dune::FVDiffusionVelocity<GridType, NumberType, VC> diffusion(grid, problem);
        //Dune::MimeticDiffusion<GridType, NumberType, VC> diffusion(grid, problem, grid.maxLevel());


        diffusion.pressure();
        std::cout << "pressure calculation took " << timer.elapsed() << " seconds" << std::endl;
        printvector(std::cout, variables.pressure, "pressure", "row", 200, 1, 3);
        variables.vtkout("fv", 0);

        diffusion.calcTotalVelocity();
        printvector(std::cout, variables.velocity, "velocity", "row", 4, 1, 3);

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
//#else
//
//int main (int argc , char **argv) try
//{
//  std::cout << "Please install the UG library." << std::endl;
//
//  return 1;
//}
//catch (...)
//{
//    std::cerr << "Generic exception!" << std::endl;
//    return 2;
//}
//#endif
