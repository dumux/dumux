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

#ifndef DUMUX_TEST2_CONSTVELPROBLEM_FRACTIONAL_FLOW_HH
#define DUMUX_TEST2_CONSTVELPROBLEM_FRACTIONAL_FLOW_HH

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/lnapl.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/porousmediumflow/2pff/implicit/propertydefaults.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include "verbosenewtoncontroller.hh"

#include "test2constvelspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class ConstVelProblem;

//////////
// Specify the properties for the ConstVel problem
//////////
namespace Properties
{
NEW_TYPE_TAG(ConstVelProblem, INHERITS_FROM(TwoPFractionalFlow, ConstVelSpatialParams));
NEW_TYPE_TAG(ConstVelCCProblem, INHERITS_FROM(CCTpfaModel, ConstVelProblem));

SET_TYPE_PROP(ConstVelCCProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ConstVelProblem, Problem, ConstVelProblem<TypeTag>);

// Set the wetting phase
SET_PROP(ConstVelProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(ConstVelProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Constant<TypeTag, Scalar> > type;
};

// Linear solver settings
SET_TYPE_PROP(ConstVelCCProblem, LinearSolver, UMFPackBackend<TypeTag> );

SET_INT_PROP(ConstVelCCProblem, Formulation, TwoPFormulation::pnsw);

SET_TYPE_PROP(ConstVelProblem, NewtonController, VerboseNewtonController<TypeTag>);

NEW_PROP_TAG(BaseProblem);
SET_TYPE_PROP(ConstVelCCProblem, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);
}

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_cc2p -parameterFile test_cc2p.input</tt>
 */
template <class TypeTag >
class ConstVelProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseProblem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using WettingPhase = typename GET_PROP_TYPE(TypeTag, WettingPhase);
    using NonwettingPhase = typename GET_PROP_TYPE(TypeTag, NonwettingPhase);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);

    enum {
        // equation indices
        transportEqIdx = Indices::transportEqIdx,
        saturationIdx = Indices::saturationIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        // world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    enum { adaptiveGrid = GET_PROP_VALUE(TypeTag, AdaptiveGrid) };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Sources = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    ConstVelProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        temperature_ = 273.15 + 20; // -> 20°C
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        using SatVector = std::array<Scalar, 2>;
        initSaturations_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, SatVector, Problem, InitialSaturations);

        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, EpisodeLength);

        this->timeManager().startNextEpisode(episodeLength_);

    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }

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
//        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos) || onLowerBoundary_(globalPos)) {
//            values.setAllNeumann();
//        }
//        else {
//            values.setAllDirichlet();
//        }
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
        PrimaryVariables values(0.0);
        values[saturationIdx] = 0.9;
        return values;
    }

    NeumannFluxes neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolvars,
                             const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);
        // const auto globalPos = scvf.ipGlobal();
        //
        // static const GlobalPosition v_t = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, Problem, TotalVelocity);
        //
        // if(onLowerBoundary_(globalPos))
        // {
        //     //we assume constant saturation at the lower boundary
        //     const auto& volVars = elemVolvars[scvf.insideScvIdx()];
        //
        //     Scalar mobW = volVars.mobility(wPhaseIdx);
        //     Scalar mobN = volVars.mobility(nPhaseIdx);
        //     Scalar mobT = mobW + mobN;
        //     auto K = this->spatialParams().permeabilityAtPos(globalPos);
        //
        //     Scalar densityW = volVars.density(wPhaseIdx);
        //     Scalar densityN = volVars.density(nPhaseIdx);
        //
        //     //ToDo use fluxVars for calculation
        //     Scalar fracW = mobW/mobT;
        //     Scalar velW = fracW*(v_t*scvf.unitOuterNormal()) + fracW*mobN*K*(densityN-densityW)*this->gravity()[dimWorld-1];
        //     values[transportEqIdx] = velW;
        // }
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
        if(globalPos[1] > 0.3*this->bBoxMax()[1] + eps_)
            values[saturationIdx] = initSaturations_[0];
        else
            values[saturationIdx] = initSaturations_[1];

        return values;
    }
    // \}

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     */
    bool shouldWriteRestartFile() const
    { return false; }

    bool shouldWriteOutput() const
    {
        //return false;
        return this->timeManager().timeStepIndex() == 0 ||
               this->timeManager().episodeWillBeFinished() ||
               this->timeManager().willBeFinished();
    }

    void episodeEnd()
    {
        this->timeManager().startNextEpisode(episodeLength_);
    }

private:

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->bBoxMax()[1] - eps_;
    }

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-7;
    std::string name_;
    std::array<Scalar, 2> initSaturations_;
    Scalar episodeLength_;
};
} //end namespace

#endif
