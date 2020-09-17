// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
/*!
 * \file
 *
 * \brief Test for throat transmissibilities
 */
#include "config.h"

#include <dune/common/float_cmp.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/transmissibility1p.hh>

namespace Dumux
{

struct MockProblem{};

struct MockElement{};

struct MockElemVolVars{};

struct MockFVElementGeometry
{
    using SubControlVolumeFace = int;
};

struct MockFluxVariablesCacheOneP
{
    using Scalar = double;
    MockFluxVariablesCacheOneP(const Throat::Shape shape)
    : shape_(shape)
    {}

    Throat::Shape shape() const
    { return shape_; }

    Scalar throatShapeFactor() const
    { return Throat::shapeFactor(shape(), throatRadius()); }

    Scalar throatCrossSectionalArea() const
    {
      if (shape_ != Throat::Shape::rectangle)
        return Throat::totalCrossSectionalArea(shape(), throatRadius());
      else
        return Throat::totalCrossSectionalAreaForRectangle(throatRadius(), 10.0*2.0*throatRadius()/*height*/);
    }

    Scalar throatLength() const
    { return 1e-2; }

    Scalar throatRadius() const
    { return 1e-3; }

    Throat::Shape throatCrossSectionShape() const
    { return shape_; }

private:
    Throat::Shape shape_;
};


int testAll(const MockFluxVariablesCacheOneP& fluxVarsCache)
{
    int success = 0;
    const auto element = MockElement{};
    const auto problem = MockProblem{};
    const auto fVElementGeometry = MockFVElementGeometry{};
    const auto elemVolVars = MockElemVolVars{};
    const auto scvf = typename MockFVElementGeometry::SubControlVolumeFace{};
    const auto shape = fluxVarsCache.shape();
    const auto resultBruus = TransmissibilityBruus<double>::singlePhaseTransmissibility(problem, element, fVElementGeometry, scvf, elemVolVars, fluxVarsCache, 0);
    const auto resultPatzekSilin = TransmissibilityPatzekSilin<double>::singlePhaseTransmissibility(problem, element, fVElementGeometry, scvf, elemVolVars, fluxVarsCache, 0);

    std::cout << "\nShape: " << Throat::shapeToString(shape) << std::endl;
    std::cout << "TransmissibilityBruus " << resultBruus << std::endl;
    std::cout << "TransmissibilityPatzekSilin " << resultPatzekSilin << std::endl;

    if (shape != Throat::Shape::square) // both laws yield slightly different values for squares
    {
        if (Dune::FloatCmp::ne<double>(resultBruus, resultPatzekSilin))
        {
            std::cout << "TransmissibilityBruus for " << Throat::shapeToString(shape) << " is wrong" <<  std::endl;
            ++success;
        }
    }
    else if (Dune::FloatCmp::ne<double>(resultBruus, resultPatzekSilin, 1e-2))
    {
        std::cout << "TransmissibilityBruus for " << Throat::shapeToString(shape) << " is wrong" <<  std::endl;
        ++success;
    }

    return success;
}

int testBruus(const MockFluxVariablesCacheOneP& fluxVarsCache, const double reference)
{
    int success = 0;
    const auto element = MockElement{};
    const auto problem = MockProblem{};
    const auto fVElementGeometry = MockFVElementGeometry{};
    const auto elemVolVars = MockElemVolVars{};
    const auto scvf = typename MockFVElementGeometry::SubControlVolumeFace{};
    const auto shape = fluxVarsCache.shape();
    const auto resultBruus = TransmissibilityBruus<double>::singlePhaseTransmissibility(problem, element, fVElementGeometry, scvf, elemVolVars, fluxVarsCache, 0);

    std::cout << "\nShape: " << Throat::shapeToString(shape) << std::endl;
    std::cout << "TransmissibilityBruus " << resultBruus << std::endl;

    if (Dune::FloatCmp::ne<double>(resultBruus, reference, 1e-5))
    {
        std::cout << "TransmissibilityBruus for " << Throat::shapeToString(shape) << " is wrong" <<  std::endl;
        ++success;
    }

    return success;
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;
    int success = 0;

    auto fluxVarsCacheCircle = MockFluxVariablesCacheOneP(Throat::Shape::circle);
    success += testAll(fluxVarsCacheCircle);

    auto fluxVarsCacheSquare= MockFluxVariablesCacheOneP(Throat::Shape::square);
    success += testAll(fluxVarsCacheSquare);

    auto fluxVarsCacheEquilateralTriangle = MockFluxVariablesCacheOneP(Throat::Shape::equilateralTriangle);
    success += testAll(fluxVarsCacheEquilateralTriangle);

    auto fluxVarsCacheTwoPlates = MockFluxVariablesCacheOneP(Throat::Shape::twoPlates);
    success += testBruus(fluxVarsCacheTwoPlates, 6.66667e-08);

    auto fluxVarsCacheRectangle = MockFluxVariablesCacheOneP(Throat::Shape::rectangle);
    success += testBruus(fluxVarsCacheRectangle, 1.24933e-09);

    if (success > 0)
        std::cout << success << " tests failed" << std::endl;

    return success;
}
