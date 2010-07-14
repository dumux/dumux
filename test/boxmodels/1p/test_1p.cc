// $Id: test_1p.cc 3779 2010-06-24 07:07:56Z bernd $
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#include "config.h"

#include "1pproblem.hh"
#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridinfo.hh>

#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

template <class Grid, class Scalar>
class CreateGrid
{
};

#if HAVE_UG
template <class Scalar>
class CreateGrid<Dune::UGGrid<2>, Scalar>
{
public:
    static Dune::UGGrid<2> *create(const Dune::FieldVector<Scalar, 2> &upperRight,
                                   const Dune::FieldVector<int, 2> &cellRes)
    {
        Dune::UGGrid<2> *grid = new Dune::UGGrid<2>;
        Dune::GridFactory<Dune::UGGrid<2> > factory(grid);
        for (int i=0; i<=cellRes[0]; i++) {
            for (int j=0; j<=cellRes[1]; j++) {
                Dune::FieldVector<double,2> pos;
                pos[0] = upperRight[0]*double(i)/cellRes[0];
                pos[1] = upperRight[1]*double(j)/cellRes[1];
                factory.insertVertex(pos);
            }
        }

        for (int i=0; i<cellRes[0]; i++) {
            for (int j=0; j<cellRes[1]; j++) {
                std::vector<unsigned int> v(4);

                int i0 = i*(cellRes[1]+1) + j;
                int i1 = i*(cellRes[1]+1) + j+1;
                int i2 = (i+1)*(cellRes[1]+1) + j;
                int i3 = (i+1)*(cellRes[1]+1) + j+1;

                v[0] = i0;
                v[1] = i1;
                v[2] = i2;
                v[3] = i3;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), v);
            }
        }

        factory.createGrid();
        grid->loadBalance();
        return grid;
    }
};

template <class Scalar>
class CreateGrid<Dune::UGGrid<3>, Scalar>
{
public:
    static Dune::UGGrid<3> *create(const Dune::FieldVector<Scalar, 3> &upperRight,
                                   const Dune::FieldVector<int, 3> &cellRes)
    {
        Dune::UGGrid<3> *grid = new Dune::UGGrid<3>;
        Dune::GridFactory<Dune::UGGrid<3> > factory(grid);
        for (int k=0; k<=cellRes[2]; k++) {
            for (int j=0; j<=cellRes[1]; j++) {
                for (int i=0; i<=cellRes[0]; i++) {
                    Dune::FieldVector<double,3> pos;
                    pos[0] = upperRight[0]*double(i)/cellRes[0];
                    pos[1] = upperRight[1]*double(j)/cellRes[1];
                    pos[2] = upperRight[2]*double(k)/cellRes[2];
                    factory.insertVertex(pos);
                }
            }
        }

        for (int k=0; k<cellRes[2]; k++) {
            for (int j=0; j<cellRes[1]; j++) {
                for (int i=0; i<cellRes[0]; i++) {
                    std::vector<unsigned int> v(8);

                    v[0] = k*(cellRes[0]+1)*(cellRes[1]+1) + j*(cellRes[0]+1) + i;
                    v[1] = k*(cellRes[0]+1)*(cellRes[1]+1) + j*(cellRes[0]+1) + i+1;
                    v[2] = k*(cellRes[0]+1)*(cellRes[1]+1) + (j+1)*(cellRes[0]+1) + i;
                    v[3] = k*(cellRes[0]+1)*(cellRes[1]+1) + (j+1)*(cellRes[0]+1) + i+1;
                    v[4] = (k+1)*(cellRes[0]+1)*(cellRes[1]+1) + j*(cellRes[0]+1) + i;
                    v[5] = (k+1)*(cellRes[0]+1)*(cellRes[1]+1) + j*(cellRes[0]+1) + i+1;
                    v[6] = (k+1)*(cellRes[0]+1)*(cellRes[1]+1) + (j+1)*(cellRes[0]+1) + i;
                    v[7] = (k+1)*(cellRes[0]+1)*(cellRes[1]+1) + (j+1)*(cellRes[0]+1) + i+1;

                    factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,3), v);
                }
            }
        }

        factory.createGrid();
        grid->loadBalance();
        return grid;
    }
};
#endif

void usage(const char *progname)
{
    std::cout << boost::format("usage: %s [--restart restartTime] gridFile.dgf tEnd dt\n")%progname;
    exit(1);
};

int main(int argc, char** argv)
{
    try {
        typedef TTAG(OnePTestProblem)                 TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar))  Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid))    Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        enum {dim = Grid::dimensionworld};
        typedef Dune::FieldVector<Scalar, dim> GlobalPosition;
        typedef Dune::GridPtr<Grid> GridPointer;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments for the program
        if (argc < 4)
            usage(argv[0]);

        // deal with the restart stuff
        int argPos = 1;
        bool restart = false;
        double restartTime = 0;
        if (std::string("--restart") == argv[argPos]) {
            restart = true;
            ++argPos;

            std::istringstream(argv[argPos++]) >> restartTime;
        }

        if (argc - argPos != 3) {
            usage(argv[0]);
        }

        double tEnd, dt;
        const char *dgfFileName = argv[argPos++];
        std::istringstream(argv[argPos++]) >> tEnd;
        std::istringstream(argv[argPos++]) >> dt;

        // create grid

        // -> load the grid from file
        GridPointer gridPtr =  GridPointer(dgfFileName);
        (*gridPtr).loadBalance();
        Dune::gridinfo(*gridPtr);

//        Grid grid;
//        Dune::GridFactory<Grid> factory(&grid);
//        Dune::GmshReader<Grid>::read(factory,"../../../dune-pdelab-howto/examples/grids/cube1045.msh",true,true);
//        factory.createGrid();
//        grid.loadBalance();

//        GlobalPosition lowerLeft(0.0);
//        GlobalPosition upperRight(1.0);
//        Dune::FieldVector<int,dim> res(8); // cell resolution
//        std::auto_ptr<Grid> grid(CreateGrid<Grid, Scalar>::create(upperRight, res));
//        Dune::gridinfo(*grid);

        // instantiate and run the concrete problem
//        Problem problem(gridPtr->levelView(gridPtr->maxLevel()));
        Problem problem(gridPtr->leafView());
        if (!problem.simulate(dt, tEnd))
            return 2;

        return 0;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
    }

    return 3;
}
