// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

#include <config.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dumux/geomechanics/stressstate/constantcementmodel.hh>
#include <dumux/io/gnuplotinterface.hh>
int main(int argc, char *argv[])
{
    using namespace Dumux;
    using CementModel = Dumux::ConstantCementModel<double>;

    // ref: http://www-odp.tamu.edu/publications/204_SR/103/103_t1.htm for SiO2
    // ref: https://surfacenet.de/calcite-217.html for CaCO3
    double phiCrit = 0.42;
    double phiB = 0.36;
    double Ks = 38e9;
    double Gs = 44e9;
    double Kc = 98e9;
    double Gc = 35e9;
    CementModel cementModel(phiCrit, phiB,
                            Ks,Gs,
                            Kc,Gc);

    std::vector<double> porosities(101);
    double phi = 0.0;
    porosities[0] = phi;
    std::generate(porosities.begin()+1,porosities.end(),[&]{return phi+=0.01;});

    std::vector<double> K(101);
    std::vector<double> G(101);

    auto fillModulus = [&](const int i, auto params){
        K[i] = params.lambda() + 2.0/3 * params.mu();
        G[i] = params.mu();
    };

    for (int i = 0; i < porosities.size(); i++) {
        auto poro = porosities[i];
        if (poro <= phiCrit)
        {
            const auto moduli = cementModel.effectiveLameModuli(poro);
            fillModulus(i,moduli);
        }
        else{
            K[i] = 0.0;
            G[i] = 0.0;
        }
    }

    GnuplotInterface<double> gnuplot{};
    gnuplot.addDataSetToPlot(porosities,K,"K");
    gnuplot.plot("K");

    gnuplot.addDataSetToPlot(porosities,G,"G");
    gnuplot.plot("G");
    return 0;
}
