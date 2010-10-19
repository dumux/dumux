#include "config.h"
#include <iostream>
#include <boost/format.hpp>
#include <iomanip>

#include <dune/grid/common/gridinfo.hh>
#include "external_interface.hh"
#include "mcwhorterproblem.hh"


int main(int argc, char** argv)
{
    try
    {
        typedef TTAG(McWhorterProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;

        static const int dim = Grid::dimension;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        //load interface-file
         Dumux::InterfaceProblemProperties interfaceProbProps("interface_MW.xml");
         double discretizationLength = interfaceProbProps.IPP_DiscretizationLength;
         int cellNumber = 2/discretizationLength;

        // define the problem dimensions
        Dune::FieldVector<Scalar, dim> lowerLeft(0);
        Dune::FieldVector<Scalar, dim> upperRight(2.0);
//        UpperRight[1] = 1;
        Dune::FieldVector<int, dim> cellNumbers(cellNumber);
//        cellNumbers[0] = 26;
        if (argc != 2)
        {
            std::cout << "usage: tEnd" << std::endl;
            return 0;
        }
        std::string arg1(argv[1]);
        std::istringstream is1(arg1);
        double tEnd;
        is1 >> tEnd;
        double dt = tEnd;
        // create a grid object
        typedef Dune::SGrid<dim, dim> GridType;
        typedef GridType::LeafGridView GridView;

        // grid reference
        GridType grid(cellNumbers,lowerLeft,upperRight);
        Dune::gridinfo(grid);


        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////

        Dumux::InterfaceFluidProperties interfaceFluidProps("interface_MW.xml");
        typedef GET_PROP_TYPE(TypeTag, PTAG(WettingPhase)) WettingPhase;
        typedef GET_PROP_TYPE(TypeTag, PTAG(NonwettingPhase)) NonwettingPhase;

        WettingPhase::Component::setViscosity(interfaceFluidProps.IFP_ViscosityWettingFluid);
        NonwettingPhase::Component::setViscosity(interfaceFluidProps.IFP_ViscosityNonWettingFluid);

        Problem problem(grid.leafView(), lowerLeft, upperRight);

        problem.timeManager().init(problem, 0, dt, tEnd);
        problem.timeManager().run();

        return 0;

    } catch (Dune::Exception &e)
    {
        std::cerr << "Dune reported error: " << e << std::endl;
    } catch (...)
    {
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
