// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <dumux/common/parameters.hh>
#include <dune/common/fvector.hh>
#include <dumux/io/container.hh>
#include "analyticalsolution.hh"

int main(int argc, char * argv[])
{
    using namespace Dumux;

    // initialize parameter tree
    Parameters::init(argc, argv);
    const auto posNormx = readFileToContainer<std::vector<double>>("pos.txt");
    MandelAnalyticalSolution<double> solver;

    std::vector<int> timeSeries{10,50, 100,1000,5000,8000,10000,20000,30000,50000};
    for(const auto & t : timeSeries)
    {
        std::vector<double> pSol;
        std::vector<double> uxSol;

        for(auto uNormx : posNormx)
        {
            auto v = Dune::FieldVector<double,2>{uNormx*100,0.0};
            auto p = solver.pressureAtPos(v,t) * 100 / 6e8;
            auto ux = solver.displacementAtPos(v,t)[0];
            pSol.push_back(p);
            uxSol.push_back(ux/100);
           // std::cout <<"uNormx"<< uNormx << "p= "<< p<< std::endl;
        }
        writeContainerToFile(pSol,"pSol"+std::to_string(t)+".txt");
        writeContainerToFile(uxSol,"uxSol"+std::to_string(t)+".txt");
    }

    return 0;
}
