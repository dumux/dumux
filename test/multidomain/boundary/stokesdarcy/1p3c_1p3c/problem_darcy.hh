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
 * \file
 * \ingroup BoundaryTests
 * \brief A Darcy test problem using Maxwell-Stefan diffusion.
 */

#ifndef DUMUX_DARCY_SUBPROBLEM_ONEPTHREEC_HH
#define DUMUX_DARCY_SUBPROBLEM_ONEPTHREEC_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/flux/maxwellstefanslaw.hh>

#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>

#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "../1p2c_1p2c/spatialparams.hh"
#include "h2n2co2fluidsystem.hh"

namespace Dumux {

template <class TypeTag>
class DarcySubProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct DarcyOnePThreeC { using InheritsFrom = std::tuple<OnePNC, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOnePThreeC> { using type = Dumux::DarcySubProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyOnePThreeC> { using type = FluidSystems::H2N2CO2FluidSystem<GetPropType<TypeTag, Properties::Scalar>>; };

// Use moles
template<class TypeTag>
struct UseMoles<TypeTag, TTag::DarcyOnePThreeC> { static constexpr bool value = true; };

// Do not replace one equation with a total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::DarcyOnePThreeC> { static constexpr int value = 3; };

//! Use a model with constant tortuosity for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::DarcyOnePThreeC>
{ using type = DiffusivityConstantTortuosity<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyOnePThreeC> { using type = Dune::YaspGrid<2>; };

//Set the diffusion type
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::DarcyOnePThreeC> { using type = MaxwellStefansLaw<TypeTag>; };

// Set the spatial paramaters type
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyOnePThreeC>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Properties

template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        // grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld,

        // primary variable indices
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx,
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    DarcySubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Darcy"), eps_(1e-7), couplingManager_(couplingManager)
    {
        pressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Pressure");
        initialMoleFractionN2_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialMoleFractionN2");
        initialMoleFractionCO2_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialMoleFractionCO2");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
      * \brief Specifies which kind of boundary condition should be
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scvf The boundary sub control volume face
      */
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values.setAllCouplingNeumann();


        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The boundary sub control volume face
     *
     * For this method, the \a values variable stores primary variables.
     */
    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVarsCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if (couplingManager().isCoupledEntity(couplingManager().darcyIdx, scvf))
            values = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, scvf, DiffusionCoefficientAveragingType::harmonic);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub control volume.
     *
     * \param element The element for which the source term is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scv The sub control volume
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    { return NumEqVector(0.0); }

    // \}

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = pressure_;
        values[conti0EqIdx + 1] = initialMoleFractionN2_;
        values[conti0EqIdx + 2] = initialMoleFractionCO2_;

        return values;
    }

    // \}

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    Scalar eps_;
    Scalar pressure_;
    Scalar initialMoleFractionN2_;
    Scalar initialMoleFractionCO2_;


    std::shared_ptr<CouplingManager> couplingManager_;
};
} // end namespace Dumux

#endif
