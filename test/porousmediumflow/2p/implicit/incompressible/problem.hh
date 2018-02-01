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
 * \ingroup TwoPTests
 * \brief The properties for the incompressible 2p test
 */
#ifndef DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROBLEM_HH
#define DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROBLEM_HH

#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>

#include <dumux/material/components/dnapl.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/2p/incompressiblelocalresidual.hh>

#include "spatialparams.hh"

namespace Dumux
{
// forward declarations
template<class TypeTag> class TwoPTestProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPIncompressible, INHERITS_FROM(TwoP));
NEW_TYPE_TAG(TwoPIncompressibleTpfa, INHERITS_FROM(CCTpfaModel, TwoPIncompressible, SpatialParams));
NEW_TYPE_TAG(TwoPIncompressibleMpfa, INHERITS_FROM(CCMpfaModel, TwoPIncompressible, SpatialParams));
NEW_TYPE_TAG(TwoPIncompressibleBox, INHERITS_FROM(BoxModel, TwoPIncompressible, SpatialParams));

// Set the grid type
SET_TYPE_PROP(TwoPIncompressible, Grid, Dune::YaspGrid<2>);

// Set the problem type
SET_TYPE_PROP(TwoPIncompressible, Problem, TwoPTestProblem<TypeTag>);

// the local residual containing the analytic derivative methods
SET_TYPE_PROP(TwoPIncompressible, LocalResidual, TwoPIncompressibleLocalResidual<TypeTag>);

// Set the fluid system
SET_PROP(TwoPIncompressible, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using WettingPhase = FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::LiquidPhase<Scalar, DNAPL<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Enable caching
SET_BOOL_PROP(TwoPIncompressible, EnableGridVolumeVariablesCache, false);
SET_BOOL_PROP(TwoPIncompressible, EnableGridFluxVariablesCache, false);
SET_BOOL_PROP(TwoPIncompressible, EnableFVGridGeometryCache, false);
} // end namespace Properties

/*!
 * \ingroup TwoPTests
 * \brief The incompressible 2p test problem.
 */
template<class TypeTag>
class TwoPTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,
        contiNEqIdx = Indices::contiNEqIdx
    };

public:
    TwoPTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry) {}

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
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
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
        fluidState.setTemperature(temperature());
        fluidState.setPressure(FluidSystem::wPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(FluidSystem::nPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, FluidSystem::wPhaseIdx);

        Scalar height = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
        Scalar depth = this->fvGridGeometry().bBoxMax()[1] - globalPos[1];
        Scalar alpha = 1 + 1.5/height;
        Scalar width = this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0];
        Scalar factor = (width*alpha + (1.0 - alpha)*globalPos[0])/width;

        // hydrostatic pressure scaled by alpha
        values[pwIdx] = 1e5 - factor*densityW*this->gravity()[1]*depth;
        values[snIdx] = 0.0;

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);
        if (onInlet_(globalPos))
            values[contiNEqIdx] = -0.04; // kg / (m * s)
        return values;
    }

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
        fluidState.setTemperature(temperature());
        fluidState.setPressure(FluidSystem::wPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(FluidSystem::nPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, FluidSystem::wPhaseIdx);

        Scalar depth = this->fvGridGeometry().bBoxMax()[1] - globalPos[1];

        // hydrostatic pressure
        values[pwIdx] = 1e5 - densityW*this->gravity()[1]*depth;
        values[snIdx] = 0.0;
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
