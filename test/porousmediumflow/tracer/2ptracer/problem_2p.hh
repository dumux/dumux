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
/*!
 * \ingroup TracerTests
 * \brief The properties for the incompressible 2p test
 */
#ifndef DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROBLEM_HH
#define DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/2p/incompressiblelocalresidual.hh>

#include "spatialparams_2p.hh"

#ifndef ENABLEINTERFACESOLVER
#define ENABLEINTERFACESOLVER 0
#endif

namespace Dumux {
// forward declarations
template<class TypeTag> class TwoPTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct TwoPIncompressible { using InheritsFrom = std::tuple<TwoP>; };
struct TwoPIncompressibleTpfa { using InheritsFrom = std::tuple<TwoPIncompressible, CCTpfaModel>; };
struct TwoPIncompressibleBox { using InheritsFrom = std::tuple<TwoPIncompressible, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPIncompressible> { using type = Dune::YaspGrid<2>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPIncompressible> { using type = TwoPTestProblem<TypeTag>; };

// the local residual containing the analytic derivative methods
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::TwoPIncompressible> { using type = TwoPIncompressibleLocalResidual<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPIncompressible>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPIncompressible>
{
private:
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = TwoPTestSpatialParams<FVGridGeometry, Scalar>;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TwoPIncompressible> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TwoPIncompressible> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::TwoPIncompressible> { static constexpr bool value = true; };

// Maybe enable the box-interface solver
template<class TypeTag>
struct EnableBoxInterfaceSolver<TypeTag, TTag::TwoPIncompressible> { static constexpr bool value = ENABLEINTERFACESOLVER; };
} // end namespace Properties

/*!
 * \ingroup TracerTests
 * \brief The incompressible 2p test problem.
 */
template<class TypeTag>
class TwoPTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        pressureH2OIdx = Indices::pressureIdx,
        saturationDNAPLIdx = Indices::saturationIdx,
        contiDNAPLEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx,
        waterPhaseIdx = FluidSystem::phase0Idx,
        dnaplPhaseIdx = FluidSystem::phase1Idx
    };

public:
    TwoPTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry) {}

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if (onLowerBoundary_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        GetPropType<TypeTag, Properties::FluidState> fluidState;
        fluidState.setTemperature(temperature());
        fluidState.setPressure(waterPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(dnaplPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, waterPhaseIdx);
        Scalar depth = this->fvGridGeometry().bBoxMax()[1] - globalPos[1];
        // hydrostatic pressure
        values[pressureH2OIdx] = 1e5 - densityW*this->gravity()[1]*depth;
        values[saturationDNAPLIdx] = 0.0;

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);
        if (onInlet_(globalPos))
        {
            values[contiDNAPLEqIdx] = -0.05; // kg / (m * s)
            values[Indices::conti0EqIdx] = -0.05;
        }

        // in the test with the oil wet lens, use higher injection rate
        if (this->spatialParams().lensIsOilWet())
            values[contiDNAPLEqIdx] *= 10;

        return values;
    }

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        GetPropType<TypeTag, Properties::FluidState> fluidState;
        fluidState.setTemperature(temperature());
        fluidState.setPressure(waterPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(dnaplPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, waterPhaseIdx);

        Scalar depth = this->fvGridGeometry().bBoxMax()[1] - globalPos[1];

        // hydrostatic pressure
        values[pressureH2OIdx] = 1e5 - densityW*this->gravity()[1]*depth;
        values[saturationDNAPLIdx] = 0;
        return values;
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0];
        Scalar lambda = (this->fvGridGeometry().bBoxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
