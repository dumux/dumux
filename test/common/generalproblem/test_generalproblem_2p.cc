// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief test for the two-phase box model
 */
#include "config.h"

#include "generallensproblem.hh"

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <iostream>


////////////////////////
// the main function
////////////////////////
void usage(const char *progname)
{
    std::cout << "usage: " << progname << " --box/--decoupled tEnd dt [refineLevel]\n";
    exit(1);
}

int main(int argc, char** argv)
{
#ifdef NDEBUG
    try {
#endif
        //TypeTag which chooses the box model
        typedef TTAG(BoxGeneralLensProblem) BoxTypeTag;
        //TypeTag which chooses the decoupled model
        typedef TTAG(DecoupledGeneralLensProblem) DecoupledTypeTag;
        typedef GET_PROP_TYPE(BoxTypeTag, PTAG(Scalar)) Scalar;
        typedef GET_PROP_TYPE(BoxTypeTag, PTAG(Grid)) Grid;
//        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef GET_PROP_TYPE(BoxTypeTag, PTAG(TimeManager)) BoxTimeManager;
        typedef GET_PROP_TYPE(DecoupledTypeTag, PTAG(TimeManager)) DecoupledTimeManager;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;

        static const int dim = Grid::dimension;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        if (argc < 4)
            usage(argv[0]);

        int argPos = 1;
        bool useBoxModel = false;
        bool useDecoupledModel = false;
        if (std::string(argv[argPos]) == "--box")
        {
            useBoxModel = true;
        }
        else if (std::string(argv[argPos]) == "--decoupled")
        {
            useDecoupledModel = true;
        }
        ++argPos;

        // read the initial time step and the end time
        double tEnd, dt;
        std::istringstream(argv[argPos++]) >> tEnd;
        std::istringstream(argv[argPos++]) >> dt;

        if (argc - argPos > 1) {
            usage(argv[0]);
        }

        double refineLevel = 0;
        if (argc - argPos == 1)
        std::istringstream(argv[argPos]) >> refineLevel;

        ////////////////////////////////////////////////////////////
        // create the grid
        ////////////////////////////////////////////////////////////
        GlobalPosition lowerLeftCorner(0.0);
        GlobalPosition domainSize;
        Dune::array< unsigned int, dim > numberOfCells; // cell resolution
        domainSize[0] = 6.0;
        domainSize[1] = 4.0;

        numberOfCells[0] = 48;
        numberOfCells[1] = 32;

        Dune::shared_ptr<Grid> grid(Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeftCorner, domainSize, numberOfCells));
        grid->globalRefine(refineLevel);

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////

        // specify dimensions of the low-permeable lens
        GlobalPosition lowerLeftLens, upperRightLens;
        lowerLeftLens[0] = 1.0;
        lowerLeftLens[1] = 2.0;
        upperRightLens[0] = 4.0;
        upperRightLens[1] = 3.0;

        // instantiate and run the concrete problem
        if (useBoxModel)
        {
            BoxTimeManager timeManager;
            Dumux::GeneralLensProblem<BoxTypeTag> problem(timeManager, grid->leafView(), lowerLeftLens, upperRightLens);
            problem.setName("generallens_box");
            timeManager.init(problem, 0, dt, tEnd, false);
            timeManager.run();
            return 0;
        }
        else if (useDecoupledModel)
        {
            DecoupledTimeManager timeManager;
            Dumux::GeneralLensProblem<DecoupledTypeTag> problem(timeManager, grid->leafView(), lowerLeftLens, upperRightLens);
            problem.setName("generallens_decoupled");
            timeManager.init(problem, 0, dt, tEnd, false);
            timeManager.run();
            return 0;
        }
        else
        {
            std::cout<<"No valid model chosen!\n";
            usage(argv[0]);
        }
#ifdef NDEBUG
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        throw;
    }
#endif // NDEBUG

    return 3;
}
