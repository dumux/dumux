// $Id:$
/*****************************************************************************
 *   Copyright (C) <YEARS> by <ADD_AUTHOR_HERE>                              *
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
// std lib includes:
#include <iostream>
#include <iomanip>
// dune stuff:
#include <dune/grid/sgrid.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
// dumux time discretization:
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/expliciteulerstep.hh"
// dumux material properties:
#include "dumux/material/fluids/water_air.hh"
#include <dumux/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>
// problem definition and model:
#include "dumux/transport/problems/testproblem_2p2c.hh"
#include "dumux/transport/fv/decoupled2p2c.hh"


/*********************************************
 *  Jochen Fritz, 2009                       *
 *  Test application to class Decoupled2p2c  *
 *********************************************/

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=3;

    // create a grid object
    typedef double Scalar;

    typedef Dune::SGrid<dim,dim> Grid;
    typedef Grid::LeafGridView GridView;
    typedef Dune::FieldVector<Grid::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(10);
    FieldVector L(0);
    FieldVector H(10);
    Grid grid(N,L,H);
    GridView gridview(grid.leafView());

    double tStart = 0;
    double tEnd = 3e4;
    int modulo = 1;
    double cFLFactor = 0.7;

    Dune::Liq_WaterAir wetmat;
    Dune::Gas_WaterAir nonwetmat;

    Dune::HomogeneousSoil<Grid, Scalar> soil;

    Dune::TwoPhaseRelations<Grid, Scalar> materialLaw(soil, wetmat, nonwetmat);

    Dune::VariableClass2p2c<GridView,Scalar> var(gridview);

    typedef Dune::Testproblem_2p2c<GridView, Scalar> TransProb;
    TransProb problem(gridview, var, wetmat, nonwetmat, soil, grid.maxLevel(), materialLaw, false);

    typedef Dune::Decoupled2p2c<GridView, Scalar> ModelType;
    ModelType model(gridview, problem);

    Dune::ExplicitEulerStep<Grid, ModelType> timestep;
    Dune::TimeLoop<Grid, ModelType > timeloop(tStart, tEnd, "2p2c", modulo, cFLFactor, 1e100, 1e100, timestep);

    Dune::Timer timer;
    timer.reset();
    timeloop.execute(model);
    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

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
