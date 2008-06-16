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
 *
 * This has been modeled to match as closely as possible to the one
 * described at:
 *
 * Sheta, Hussam: "Simulation von Mehrphasenvorgaengen in poroesen
 *     Medien unter Einbeziehung von Hystereseeffekten", PhD theses,
 *     Braunschweig 1999, pp. 112
 */

#include "StupidConfig.hh"
#include "LenhardConfig.hh"

#include "LenhardProblem.hh"

#include <iostream>
#include <boost/format.hpp>

int main(int argc, char** argv)
{
    try{
        typedef double Scalar;

        if (argc != 3) {
            std::cout << boost::format("usage: %s tEnd dt\n")%argv[0];
            return 1;
        }
    double tEnd, dt;
    std::istringstream(argv[1]) >> tEnd;
    std::istringstream(argv[2]) >> dt;

        Stupid::Lenhard::PwSnLenhardProblem<Scalar> problem(dt, tEnd);
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
