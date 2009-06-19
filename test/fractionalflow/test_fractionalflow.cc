// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch, Jochen Fritz, Markus Wolff   *
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
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "test_fractionalflow_soilproperties.hh"
#include <dumux/material/fluids/water.hh>
#include <dumux/material/fluids/oil.hh>
#include <dumux/material/twophaserelations.hh>
#include "test_fractionalflow_testproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/diffusion/fv/fvtotalvelocity2p.hh"
#include "dumux/transport/fv/fvsaturationwetting2p.hh"
#include "dumux/transport/fv/capillarydiffusion.hh"
#include "dumux/fractionalflow/impes/impes.hh"
#include <dune/disc/operators/boundaryconditions.hh>
#include "dumux/timedisc/expliciteulerstep.hh"
#include "dumux/fractionalflow/variableclass2p.hh"


int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=2;

        if (argc != 2) {
            std::cout << "\nUsage: test_fractionalflow entryPressure\n" << std::endl;
            std::cout << "- entryPressure: the capillary entry pressure >= 0." << std::endl;
            std::cout << "  If 0, no capillary effects are taken into account.\n" << std::endl;
            return (1);
        }

        double entryPressure = 0;
        std::string arg1(argv[1]);
        std::istringstream is1(arg1);
        is1 >> entryPressure;

        typedef double Scalar;
        typedef Dune::SGrid<dim,dim> Grid;
        typedef Grid::LevelGridView GridView;
        typedef Dune::FieldVector<Grid::ctype,dim> FieldVector;

        Dune::FieldVector<int,dim> N(1); N[0] = 64;
        FieldVector L(0);
        FieldVector H(60); H[0] = 300;
        Grid grid(N,L,H);

        grid.globalRefine(0);
        GridView gridView(grid.levelView(0));

        Dune::Water wetmat;
        Dune::Oil nonwetmat;

        Dune::FractionalFlowTestSoil<Grid, Scalar> soil(entryPressure);

        Dune::TwoPhaseRelations<Grid, Scalar> materialLaw(soil, wetmat, nonwetmat);

        typedef Dune::VariableClass<GridView, Scalar> VariableType;

        VariableType variables(gridView);

        typedef Dune::FractionalFlowTestProblem<GridView, Scalar, VariableType> Problem;
        Problem problem(variables, wetmat, nonwetmat, soil, materialLaw,L, H);

        typedef Dune::FVTotalVelocity2P<GridView, Scalar, VariableType, Problem> DiffusionType;
        DiffusionType diffusion(gridView, problem, "pglobal");

        Dune::CapillaryDiffusion<GridView, Scalar, VariableType, Problem> capillaryDiffusion(problem, soil);

        typedef Dune::FVSaturationWetting2P<GridView, Scalar, VariableType, Problem> TransportType;
        TransportType transport(gridView, problem, "vt", capillaryDiffusion);

        int iterFlag = 0;
        int nIter = 30;
        double maxDefect = 1e-5;
        typedef Dune::IMPES<GridView, DiffusionType, TransportType, VariableType> IMPESType;
        IMPESType impes(diffusion, transport, iterFlag, nIter, maxDefect);

        double tStart = 0;
        double tEnd = 5e7;
        const char* fileName = "test_fractionalflow";
        int modulo = 500;
        double cFLFactor = 0.5;
        Dune::TimeLoop<Grid, IMPESType > timeloop(tStart, tEnd, fileName, modulo, cFLFactor);

        timeloop.execute(impes);

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
    }
}
