// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \ingroup TwoPTests
 * \brief The properties for the incompressible 2p test.
 */
#ifndef DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROBLEM_HH
#define DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief The incompressible 2p test problem.
 */
template<class TypeTag>
class TwoPTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        pressureH2OIdx = Indices::pressureIdx,
        saturationDNAPLIdx = Indices::saturationIdx,
        contiDNAPLEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx,
        waterPhaseIdx = FluidSystem::phase0Idx,
        dnaplPhaseIdx = FluidSystem::phase1Idx
    };

public:
    TwoPTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry) {}

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
        fluidState.setTemperature(this->spatialParams().temperatureAtPos(globalPos));
        fluidState.setPressure(waterPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(dnaplPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, waterPhaseIdx);

        Scalar height = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];
        Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        Scalar alpha = 1 + 1.5/height;
        Scalar width = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
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
        fluidState.setTemperature(this->spatialParams().temperatureAtPos(globalPos));
        fluidState.setPressure(waterPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(dnaplPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, waterPhaseIdx);

        Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];

        // hydrostatic pressure
        values[pressureH2OIdx] = 1e5 - densityW*this->spatialParams().gravity(globalPos)[1]*depth;
        values[saturationDNAPLIdx] = 0;
        return values;
    }


private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        Scalar lambda = (this->gridGeometry().bBoxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
