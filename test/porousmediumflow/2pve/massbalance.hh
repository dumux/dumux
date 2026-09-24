// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVETests
 * \brief Helper functions to evaluate the mass balance of the gas phase in the domain. The evaluation is only valid while the gas plume tip does not reach the right boundary.
 */

#ifndef DUMUX_TEST_TWOPVE_MASSBALANCE_HH
#define DUMUX_TEST_TWOPVE_MASSBALANCE_HH

#include <algorithm>
#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/scvandscvfiterators.hh>

namespace Dumux::VETest {

template<typename Scalar>
struct VEMassBalance
{
    Scalar nonwettingMassCoarse{};
    Scalar nonwettingMassFine{};
    Scalar expectedInjectedMass{};

    Scalar relativeErrorCoarse() const
    { return (expectedInjectedMass - nonwettingMassCoarse)/expectedInjectedMass; }

    Scalar relativeErrorFine() const
    { return (expectedInjectedMass - nonwettingMassFine)/expectedInjectedMass; }
};

/*!
 * \brief Computes the gas-phase mass on the coarse and fine level, as well as the expected gas mass in the domain
 *
 * \param fvGridGeometryVE     coarse-level grid geometry
 * \param fvGridGeometryVEFine fine-level grid geometry
 * \param solution             coarse-level solution vector
 * \param problemVE            coarse-level problem
 * \param timeLoop             container with time properties
 */
template<typename TypeTag>
auto computeMassBalance(const GetPropType<TypeTag, Properties::GridGeometry>& fvGridGeometryVE,
                        const GetPropType<TypeTag, Properties::GridGeometry>& fvGridGeometryVEFine,
                        const GetPropType<TypeTag, Properties::SolutionVector>& solution,
                        const GetPropType<TypeTag, Properties::Problem>& problemVE,
                        const TimeLoop<GetPropType<TypeTag, Properties::Scalar>>& timeLoop)
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NonwettingPhase = typename GetProp<TypeTag, Properties::FluidSystem>::NonwettingPhase;

    constexpr int dim = GridView::dimension;
    auto injectionRate = getParam<Scalar>("BoundaryConditions.InjectionRate");
    const GlobalPosition lowerLeft = getParam<GlobalPosition>("Grid.LowerLeft");
    const GlobalPosition upperRight = getParam<GlobalPosition>("Grid.UpperRight");
    const auto dummyElement = *(fvGridGeometryVE.gridView().template begin<0>());
    const Scalar dummyTemperature = problemVE.spatialParams().temperatureAtPos(dummyElement.geometry().center());
    const auto& fineSol = problemVE.getFineLevelView()->solution();
    const auto& spatialParamsCoarse = problemVE.spatialParams();
    const auto& spatialParamsFine = problemVE.getFineLevelView()->spatialParams();
    Scalar domainHeight = upperRight[dim-1] - lowerLeft[dim-1];
    auto massBalance = VEMassBalance<Scalar>();

    // iteration over coarse-level elements
    for (const auto& element : Dune::elements(fvGridGeometryVE.gridView()))
    {
        auto elementIdx = fvGridGeometryVE.elementMapper().index(element);
        auto satNwCoarse = solution[elementIdx][1];
        auto pwCoarse = solution[elementIdx][0];
        auto densityNw = NonwettingPhase::density(dummyTemperature, pwCoarse);
        auto fvGeometryVE = localView(fvGridGeometryVE);
        fvGeometryVE.bind(element);
        Scalar coarsePorosity = spatialParamsCoarse.porosityCoarseAtElement(element);

        for (const auto& scvVE : scvs(fvGeometryVE))
        {
            Scalar realElementVolumeCoarse = scvVE.volume()/domainHeight;
            massBalance.nonwettingMassCoarse += coarsePorosity * densityNw * satNwCoarse * realElementVolumeCoarse;
        }

        // iteration over fine-level elements
        auto column = problemVE.getFineLevelView()->columnMap().column(elementIdx);
        for(const auto& fineElement : column)
        {
            auto elementIdxFine = fvGridGeometryVEFine.elementMapper().index(fineElement);
            auto densityNwFine = NonwettingPhase::density(dummyTemperature, pwCoarse);
            auto satNwFine = fineSol[elementIdxFine][1];
            auto fvGeometryVEFine = localView(fvGridGeometryVEFine);
            fvGeometryVEFine.bind(fineElement);
            Scalar finePorosity = spatialParamsFine.porosityAtElement(fineElement);

            for (const auto& scv : scvs(fvGeometryVEFine))
            {
                massBalance.nonwettingMassFine += finePorosity * densityNwFine * satNwFine * scv.volume();
            }
        }
    }

    const Scalar injectionDuration = std::min(timeLoop.time(), problemVE.injectionEndTime());
    massBalance.expectedInjectedMass = injectionDuration * domainHeight * (-injectionRate);

    return massBalance;
}

/*!
 * \brief Print the state of the mass balance object
 *
 * \param massBalance object that contains current gas mass in system
 */
template<typename Scalar>
void printMassBalance(const VEMassBalance<Scalar>& massBalance)
{
    std::cout << "-------info about mass conservation-------" << std::endl;
    std::cout << "Expected injected gas mass is: " << massBalance.expectedInjectedMass << "." << std::endl;
    std::cout << "Gas mass in VE coarse system is: " << massBalance.nonwettingMassCoarse <<  " and in VE fine system: " << massBalance.nonwettingMassFine << "." << std::endl;
    std::cout << "Error on coarse VE level: " << massBalance.relativeErrorCoarse() << ", error on fine VE level: " << massBalance.relativeErrorFine() << std::endl;
    std::cout << "-----------------info end-----------------\n";
}

/*!
 * \brief Throws if the gas mass on the coarse or on the fine level deviates from the injected gas mass by more than the given relative tolerance
 *
 * \param massBalance       object that contains current gas mass in system
 * \param relativeTolerance maximum relative deviation from the injected gas mass
 */
template<typename Scalar>
void checkMassBalance(const VEMassBalance<Scalar>& massBalance, const Scalar relativeTolerance)
{
    using std::abs;
    if (abs(massBalance.relativeErrorCoarse()) > relativeTolerance || abs(massBalance.relativeErrorFine()) > relativeTolerance)
        DUNE_THROW(Dune::Exception, "Gas mass deviates from the injected gas mass by more than " << relativeTolerance
                                    << " (relative error on coarse level: " << massBalance.relativeErrorCoarse()
                                    << ", on fine level: " << massBalance.relativeErrorFine() << ")");
}

} // end namespace Dumux::VETest

#endif
