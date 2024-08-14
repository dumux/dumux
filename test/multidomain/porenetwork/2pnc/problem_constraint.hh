// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief A test problem for setting a throat constraint for the two-phase pore network model.
 */
#ifndef DUMUX_PNM2PNC_CONSTRAINT_PROBLEM_HH
#define DUMUX_PNM2PNC_CONSTRAINT_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porenetwork/2p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

template<class TypeTag>
class PNMConstraintProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename GridGeometry::LocalView::Element::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;

public:
    PNMConstraintProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                         const std::string& paramGroup = "",
                         std::shared_ptr<CouplingManager> couplingManager = nullptr)
    : ParentType(gridGeometry, paramGroup),
      couplingManager_(couplingManager)
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return { 0.0 }; }

    PrimaryVariables initial(const Element& element) const
    { return {(Scalar) couplingManager_->problem(CouplingManager::poreNetworkIndex).initialInvasionState(element)}; }

    Scalar volumeConstraint(const Element &element,
                            const FVElementGeometry &fvGeometry,
                            const ElementVolumeVariables &elemVolVars,
                            const SubControlVolume &scv) const
    {
        auto theta = elemVolVars[scv].theta();
        const auto& elemVolVarsPNM = couplingManager_->elemVolVars(element);
        using std::max; using std::min;

        const auto& state = couplingManager_->gridVariables(CouplingManager::poreNetworkIndex).gridFluxVarsCache().invasionState();
        const auto prevInvaded = state.invaded(element);

        if(!prevInvaded)
        {
            auto dp = max(elemVolVarsPNM[0].capillaryPressure(),
                          elemVolVarsPNM[1].capillaryPressure()) / couplingManager_->pcEntry(element) - 1.0;

            return ((1-theta)*max(0.0,dp) - (theta)*min(0.0,dp));
        }
        else
        {
            auto pcSnapoff = couplingManager_->pcSnapoff(element);
            auto dp = min(elemVolVarsPNM[0].capillaryPressure(),
                          elemVolVarsPNM[1].capillaryPressure()) / abs(pcSnapoff) - sign(pcSnapoff);

            return ((1-theta)*max(0.0,dp) - (theta)*min(0.0,dp));
        }
    }

private:
    std::shared_ptr<CouplingManager> couplingManager_;

};
} // end namespace Dumux

#endif
