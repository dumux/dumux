/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief Solves the "lenhard problem" using a fully coupled Pw-Sn formulation.
 */

/////////////////////////////////
// Coarse configuration parameters for the lenhard problem.
/////////////////////////////////

// if defined, most material parameters (capillary pressure, permeability,
// etc) are defined on the vertices of FE grid.
#define USE_NODE_PARAMETERS

// if defined, interface conditions of capillary pressure und and
// relative permeability are used at medium interfaces for element based
// material parameters.
//#define USE_INTERFACE_CONDITION

// use parker-lenhard hysteresis for the simulation. if undefined only
// van-genuchten without hysteresis is used.
#define USE_HYSTERESIS

// write the convergence behaviour to disk. beware, this results in a
// HUGE number of files.
//#define LENHARD_WRITE_NEWTON_STEPS

// simulate experiment 1, or experiment 2?
#define LENHARD_EXPERIMENT 2
/////////////////////////////////
// End of coarse configuration parameters for the lenhard problem.
/////////////////////////////////

#include "config.h"
#include "lenhardproblem.hh"

#include <iostream>
#include <boost/format.hpp>

int main(int argc, char** argv)
{
    try{
        typedef double Scalar;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        if (argc != 3) {
            std::cout << boost::format("usage: %s tEnd dt\n")%argv[0];
            return 1;
        }
        Scalar tEnd, dt;
        std::istringstream(argv[1]) >> tEnd;
        std::istringstream(argv[2]) >> dt;

        Dune::Lenhard::PwSnLenhardProblem<Scalar> problem(dt, tEnd);
        if (!problem.simulate())
            return 2;
        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!\n";
    }
    return 3;
}
