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
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_CONSTVELPROBLEM_FRACTIONAL_FLOW_HH
#define DUMUX_CONSTVELPROBLEM_FRACTIONAL_FLOW_HH

#include <dumux/common/parameters.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/porousmediumflow/2pff/implicit/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include "constvelspatialparams.hh"

namespace Dumux {

template <class TypeTag>
class ConstVelProblem;

//////////
// Specify the properties for the ConstVel problem
//////////
namespace Properties {
NEW_TYPE_TAG(ConstVelProblem, INHERITS_FROM(TwoPFractionalFlow));
NEW_TYPE_TAG(ConstVelCCProblem, INHERITS_FROM(CCTpfaModel, ConstVelProblem));

SET_TYPE_PROP(ConstVelProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ConstVelProblem, Problem, ConstVelProblem<TypeTag>);

// set spatial params
SET_TYPE_PROP(ConstVelProblem, SpatialParams, ConstVelSpatialParams<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                                                                    typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Set the wetting phase
SET_PROP(ConstVelProblem, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using LNapl = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
    using type = FluidSystems::TwoPImmiscible<Scalar, FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar>>, LNapl>;
};

SET_PROP(ConstVelCCProblem, Formulation)
{ static constexpr TwoPFormulation value = TwoPFormulation::p1s0; };

// SET_TYPE_PROP(ConstVelProblem, NewtonController, VerboseNewtonController<TypeTag>);
} // end namespace Properties

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_cc2p -parameterFile test_cc2p.input</tt>
 */
template <class TypeTag>
class ConstVelProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    enum {
        // equation indices
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        transportEqIdx = 0,
        totalvelocityEqIdx = 1,

        // phase indices
        wPhaseIdx = FluidSystem::phase0Idx,
        nPhaseIdx = FluidSystem::phase1Idx,

        // world dimension
        dimWorld = GridView::dimensionworld
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    /*!
     * \brief The constructor
     */
    ConstVelProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        temperature_ = 273.15 + 20; // -> 20°C
        initSaturations_ = getParam<std::array<Scalar, 2>>("Problem.InitialSaturations");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * This problem assumes a uniform temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

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
        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos) || onLowerBoundary_(globalPos)) {
            values.setAllNeumann();
        }
        else
        {
            values.setAllDirichlet();
        }
//        values.setAllNeumann();
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
        PrimaryVariables values(0.0);
        values[saturationIdx] = 0.0;
        values[pressureIdx] = 1.0e5;
        return values;
    }

    template<class ElementVolumeVariables>
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolvars,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);

         const auto globalPos = scvf.ipGlobal();

         static const GlobalPosition v_t = getParam<GlobalPosition>("Problem.TotalVelocity");
         if(onLowerBoundary_(globalPos))
         {
             //we assume constant saturation at the lower boundary
             const auto& volVars = elemVolvars[scvf.insideScvIdx()];

             Scalar mobW = volVars.mobility(wPhaseIdx);
             Scalar mobN = volVars.mobility(nPhaseIdx);
             Scalar mobT = mobW + mobN;
             auto K = this->spatialParams().permeabilityAtPos(globalPos);

             Scalar densityW = volVars.density(wPhaseIdx);
             Scalar densityN = volVars.density(nPhaseIdx);

             //ToDo use fluxVars for calculation
             Scalar fracW = mobW/mobT;
             Scalar fracN = mobN/mobT;
             Scalar fluxW = fracW*(v_t*scvf.unitOuterNormal()) + fracW*mobN*K*(densityN-densityW)*this->gravity()[dimWorld-1];
             Scalar fluxN = fracN*(v_t*scvf.unitOuterNormal()) - fracN*mobW*K*(densityN-densityW)*this->gravity()[dimWorld-1];
             fluxW *= scvf.area();
             fluxN *= scvf.area();

             values[transportEqIdx] = fluxW;
             values[totalvelocityEqIdx] = fluxW + fluxN;
         }

        return values;
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{


    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = 1.0e5;
        if(globalPos[1] > 0.65*this->fvGridGeometry().bBoxMax()[1] + eps_)
            values[saturationIdx] = initSaturations_[0];
        else
            values[saturationIdx] = initSaturations_[1];

        return values;
    }
    // \}
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

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-7;
    std::array<Scalar, 2> initSaturations_;
};

} // end namespace Dumux

#endif
