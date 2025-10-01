// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>

#include <dumux/io/container.hh>
#include <dumux/io/format.hh>

#include <dune/common/fvector.hh>

#include "analyticalsolution.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    MandelAnalyticalSolution<double> solver;

    const auto posx = readFileToContainer<std::vector<double>>("pos.txt");
    std::vector<int> times{ 10, 50, 100, 1000, 5000, 8000, 10000, 20000, 30000, 50000 };
    for(const auto& t : times)
    {
        std::vector<double> pSol;
        std::vector<double> uxSol;

        for (const auto x : posx)
        {
            auto pos = Dune::FieldVector<double, 2>{ x*100.0, 0.0 };
            auto p = solver.pressure(pos, t) * 100.0 / 6e8;
            auto ux = solver.displacement(pos, t)[0];
            pSol.push_back(p);
            uxSol.push_back(ux / 100.0);
        }

        writeContainerToFile(pSol, Fmt::format("p-{}.txt", t));
        writeContainerToFile(uxSol, Fmt::format("ux-{}.txt", t));
    }

    return 0;
}
