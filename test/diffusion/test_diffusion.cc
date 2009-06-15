#include "config.h"
#include <iostream>
#include <iomanip>
//#ifdef HAVE_UG
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/sgrid.hh>
//#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
//#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include <dumux/material/fluids/uniform.hh>
#include <dumux/material/matrixproperties.hh>
//
#include "dumux/diffusion/fv/fvwettingvelocity2p.hh"
//#include "dumux/diffusion/fe/fepressure2p.hh"
//
#include "dumux/diffusion/mimetic/mimeticpressure2p.hh"
//
#include "dumux/diffusion/problems/uniformproblem.hh"
#include "dumux/fractionalflow/variableclass2p.hh"

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
        typedef GridType::LevelGridView GridView;
        GridView gridView(grid.levelView(0));


        //Uniform mat;
        Dune::Uniform mat;

        Dune::HomogeneousSoil<GridType, NumberType> soil;
//        Dune::HeterogeneousSoil<GridType, NumberType> soil(grid, "permeab.dat", true);
//        printvector(std::cout, *(soil.permeability), "permeability", "row", 200, 1);
//        soil.permeability.vtkout("permeability", grid);

        Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, mat, mat);

        typedef Dune::VariableClass<GridView, NumberType> VC;

        double initsat = 0.8;

//        VC variables(gridView,initsat);
        //for fe discretisation -> pressure on the nodes!
        VC variables(gridView, dim, initsat);

        Dune::UniformProblem<GridView, NumberType, VC> problem(variables, mat, mat, soil, materialLaw, R);

        Dune::Timer timer;
        timer.reset();
//        Dune::LeafFEPressure2P<GridView, NumberType, VC> diffusion(gridView, problem);
        Dune::FVWettingPhaseVelocity2P<GridView, NumberType, VC> diffusion(gridView, problem, "pw","Sw");
//        Dune::MimeticPressure2P<GridView, NumberType, VC> diffusion(gridView, problem);


        diffusion.pressure();
        std::cout << "pressure calculation took " << timer.elapsed() << " seconds" << std::endl;
        printvector(std::cout, variables.pressure(), "pressure", "row", 200, 1, 3);
        variables.vtkout("fv", 0);

        diffusion.calculateVelocity();
        printvector(std::cout, variables.velocity(), "velocity", "row", 4, 1, 3);

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}

