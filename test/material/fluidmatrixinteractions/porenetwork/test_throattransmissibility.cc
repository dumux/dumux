// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test for throat transmissibilities
 */
#include "config.h"

#include <dune/common/float_cmp.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>

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
    MockFluxVariablesCacheOneP(const PoreNetwork::Throat::Shape shape)
    : shape_(shape)
    {}

    PoreNetwork::Throat::Shape shape() const
    { return shape_; }

    Scalar throatShapeFactor() const
    { return PoreNetwork::Throat::shapeFactor(shape(), throatInscribedRadius()); }

    Scalar throatCrossSectionalArea() const
    {
      if (shape_ != PoreNetwork::Throat::Shape::rectangle)
        return PoreNetwork::Throat::totalCrossSectionalArea(shape(), throatInscribedRadius());
      else
        return PoreNetwork::Throat::totalCrossSectionalAreaForRectangle(throatInscribedRadius(), 10.0*2.0*throatInscribedRadius()/*height*/);
    }

    Scalar throatLength() const
    { return 1e-2; }

    Scalar throatInscribedRadius() const
    { return 1e-3; }

    PoreNetwork::Throat::Shape throatCrossSectionShape() const
    { return shape_; }

private:
    PoreNetwork::Throat::Shape shape_;
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
    const auto resultBruus = PoreNetwork::TransmissibilityBruus<double>::singlePhaseTransmissibility(problem, element, fVElementGeometry, scvf, elemVolVars, fluxVarsCache, 0);
    const auto resultPatzekSilin = PoreNetwork::TransmissibilityPatzekSilin<double>::singlePhaseTransmissibility(problem, element, fVElementGeometry, scvf, elemVolVars, fluxVarsCache, 0);

    std::cout << "\nShape: " << PoreNetwork::Throat::shapeToString(shape) << std::endl;
    std::cout << "TransmissibilityBruus " << resultBruus << std::endl;
    std::cout << "TransmissibilityPatzekSilin " << resultPatzekSilin << std::endl;

    if (shape != PoreNetwork::Throat::Shape::square) // both laws yield slightly different values for squares
    {
        if (Dune::FloatCmp::ne<double>(resultBruus, resultPatzekSilin))
        {
            std::cout << "TransmissibilityBruus for " << PoreNetwork::Throat::shapeToString(shape) << " is wrong" <<  std::endl;
            ++success;
        }
    }
    else if (Dune::FloatCmp::ne<double>(resultBruus, resultPatzekSilin, 1e-2))
    {
        std::cout << "TransmissibilityBruus for " << PoreNetwork::Throat::shapeToString(shape) << " is wrong" <<  std::endl;
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
    const auto resultBruus = PoreNetwork::TransmissibilityBruus<double>::singlePhaseTransmissibility(problem, element, fVElementGeometry, scvf, elemVolVars, fluxVarsCache, 0);

    std::cout << "\nShape: " << PoreNetwork::Throat::shapeToString(shape) << std::endl;
    std::cout << "TransmissibilityBruus " << resultBruus << std::endl;

    if (Dune::FloatCmp::ne<double>(resultBruus, reference, 1e-5))
    {
        std::cout << "TransmissibilityBruus for " << PoreNetwork::Throat::shapeToString(shape) << " is wrong" <<  std::endl;
        ++success;
    }

    return success;
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;
    int success = 0;

    auto fluxVarsCacheCircle = MockFluxVariablesCacheOneP(PoreNetwork::Throat::Shape::circle);
    success += testAll(fluxVarsCacheCircle);

    auto fluxVarsCacheSquare= MockFluxVariablesCacheOneP(PoreNetwork::Throat::Shape::square);
    success += testAll(fluxVarsCacheSquare);

    auto fluxVarsCacheEquilateralTriangle = MockFluxVariablesCacheOneP(PoreNetwork::Throat::Shape::equilateralTriangle);
    success += testAll(fluxVarsCacheEquilateralTriangle);

    auto fluxVarsCacheTwoPlates = MockFluxVariablesCacheOneP(PoreNetwork::Throat::Shape::twoPlates);
    success += testBruus(fluxVarsCacheTwoPlates, 6.66667e-08);

    auto fluxVarsCacheRectangle = MockFluxVariablesCacheOneP(PoreNetwork::Throat::Shape::rectangle);
    success += testBruus(fluxVarsCacheRectangle, 1.24933e-09);

    if (success > 0)
        std::cout << success << " tests failed" << std::endl;

    return success;
}
