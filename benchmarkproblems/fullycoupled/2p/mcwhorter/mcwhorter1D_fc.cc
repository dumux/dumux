# include "config.h"

// disable the restart functionality since
// this doesn't work for 1D
#define DUMUX_NO_RESTART 1

#include <iostream>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "../problemdefinitions/mcwhorterproblem.hh"
#include "dumux/twophase/problems/uniformtwophaseproblem.hh"
#include "dumux/twophase/fv/boxpnsw_deprecated.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/brookscoreylaw_deprecated.hh"
#include "dumux/material/vangenuchtenlaw_deprecated.hh"
#include "dumux/io/vtkmultiwriter.hh"
#include"../problemdefinitions/mcwhorteranalytical.hh"

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=1;
        enum {BrooksCorey = 0, VanGenuchten = 1};
        typedef double NumberType;
        typedef GridType::ctype ctype;
        Dune::FieldVector<NumberType, dim> Left(0);
        Dune::FieldVector<NumberType, dim> Right(2.6);
        if (argc < 3) {
            std::cout << "usage: test_twophase tEnd dt" << std::endl;
            return 0;
        }
        std::string arg1(argv[1]);
        std::istringstream is1(arg1);
        int tEnd;
        is1 >> tEnd;
        std::string arg2(argv[2]);
        std::istringstream is2(arg2);
        int dt;
        is2 >> dt;

        // create a grid object
        typedef Dune::OneDGrid GridType;

        const int numberofelements = 52;
        double strfactor = 0;

        //vector with coordinates
        std::vector<ctype> coord;
        coord.resize(numberofelements+1);
        coord[0]=0;
        coord[1]=1;
        //generate coordinates for a stretched grid
        for (int i=2;i<numberofelements+1;i++){
            coord[i]=coord[i-1]+(coord[i-1]-coord[i-2])*(1+strfactor);
        }

        //scale coordinates to geometry
        for (int i=0;i<numberofelements+1;i++){
            coord[i]*=Right[0]/coord[numberofelements];
            std::cout << "coordinates =  " << coord[i] << std::endl;
        }

        const std::vector<ctype>& coordinates(coord);

        // grid reference
        GridType grid(coordinates);

        Dune::gridinfo(grid);

        Oil oil(0);
        Water water(0);
        //Dune::DeprecatedLinearLaw law(water,oil);
        Dune::DeprecatedBrooksCoreyLaw law(water, oil,2,5000);
        //Dune::DeprecatedVanGenuchtenLaw law(water, oil,3.1257,1.74e-4);

        //calculate with analytical solution
        Dune::McWWithAnalytical<GridType, NumberType> problem(grid,law, Left, Right);

        //calculate without analytical solution
        //      Dune::McWhorterProblem<GridType, NumberType> problem(law, Left, Right,/*VanGenuchten*/BrooksCorey);

        typedef Dune::DeprecatedBoxPnSw<GridType, NumberType> TwoPhase;
        TwoPhase twoPhase(grid, problem);

        Dune::TimeLoop<GridType, TwoPhase> timeloop(0, tEnd, dt,
                                                    "mcwhorter1D", 1);

        Dune::Timer timer;
        timer.reset();
        timeloop.execute(twoPhase);
        std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

        //printvector(std::cout, *twoPhase.u, "u", "row", 2, 1, 3);

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
