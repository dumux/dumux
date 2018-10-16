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
#ifndef DUMUX_IMPES_TRANSPORT_TEST_PROBLEM_HH
#define DUMUX_IMPES_TRANSPORT_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>

#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>

#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2p/sequential/saturation/model.hh>

#include "transportspatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class TwoPTransport;

namespace Properties {
NEW_TYPE_TAG(TwoPTransport, INHERITS_FROM(CCTpfaModel, Transport));
// Set the grid type
SET_TYPE_PROP(TwoPTransport, Grid, Dune::YaspGrid<2>);

// Set the problem type
SET_TYPE_PROP(TwoPTransport, Problem, TwoPTransport<TypeTag>);

// Set the fluid system
SET_PROP(TwoPTransport, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
#if PROBLEM == 2
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
#else
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
#endif
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Set the spatial parameters
SET_PROP(TwoPTransport, SpatialParams)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = TwoPTransportSpatialParams<FVGridGeometry, Scalar>;
};

// Enable caching
SET_BOOL_PROP(TwoPTransport, EnableGridVolumeVariablesCache, true);
SET_BOOL_PROP(TwoPTransport, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(TwoPTransport, EnableFVGridGeometryCache, true);
} // end namespace Properties

/*!
 * \ingroup TwoPTests
 * \brief The incompressible 2p test problem.
 */
template<class TypeTag>
class TwoPTransport : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    enum {
        transportEqIdx = Indices::transportEqIdx,
        saturationIdx = Indices::saturationIdx
    };
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    TwoPTransport(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        Scalar inletWidth = getParam<Scalar>("Problem.InletWidth", 1.0);
        GlobalPosition inletCenter = this->fvGridGeometry().bBoxMax();
        inletCenter[0] *= 0.5;

        inletLeftCoord_ = inletCenter;
        inletLeftCoord_[0] -=0.5*inletWidth;
        inletRightCoord_ = inletCenter;
        inletRightCoord_[0] +=0.5*inletWidth;

        inFlux_ = getParam<Scalar>("Problem.InjectionFlux", 1e-4);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
#if PROBLEM == 2
        BoundaryTypes values;
        if (onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos))
        {
            values.setAllNeumann();
        }
        else
        {
            values.setAllDirichlet();
        }
        return values;
#else
        BoundaryTypes values;
        if (onLeftBoundary_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
#endif
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
        values[saturationIdx] = 1.0;

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
        values[saturationIdx] = 1.0;

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

    bool isInlet(const GlobalPosition& globalPos) const
    {
        if (!onUpperBoundary_(globalPos))
            return false;

        for (int i = 0; i < dimWorld; i++)
        {
            if (globalPos[i] < inletLeftCoord_[i] - eps_)
                return false;
            if (globalPos[i] > inletRightCoord_[i] + eps_)
                return false;
        }
        return true;
    }

    static constexpr Scalar eps_ = 1e-6;
    Scalar inFlux_;
    GlobalPosition inletLeftCoord_;
    GlobalPosition inletRightCoord_;
};

} // end namespace Dumux

#endif
