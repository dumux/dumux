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
 * \ingroup TwoPTests
 * \brief The properties for the incompressible 2p test.
 */
#ifndef DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROBLEM_HH
#define DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROBLEM_HH

#include <dumux/porousmediumflow/problem.hh>

#include "spatialparams.hh"

namespace Dumux {

/*!
 * \ingroup TwoPTests
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
        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
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

        Scalar height = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
        Scalar depth = this->fvGridGeometry().bBoxMax()[1] - globalPos[1];
        Scalar alpha = 1 + 1.5/height;
        Scalar width = this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0];
        Scalar factor = (width*alpha + (1.0 - alpha)*globalPos[0])/width;

        // hydrostatic pressure scaled by alpha
        values[pressureH2OIdx] = 1e5 - factor*densityW*this->spatialParams().gravity(globalPos)[1]*depth;
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
            values[contiDNAPLEqIdx] = -0.04; // kg / (m * s)

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
        values[pressureH2OIdx] = 1e5 - densityW*this->spatialParams().gravity(globalPos)[1]*depth;
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
        return 293.15; // 10°C
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
