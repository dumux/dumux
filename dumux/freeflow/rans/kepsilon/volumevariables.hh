// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::KEpsilonMassVolumeVariables
 */
#ifndef DUMUX_RANS_KEPSILON_MASS_VOLUME_VARIABLES_HH
#define DUMUX_RANS_KEPSILON_MASS_VOLUME_VARIABLES_HH

#include <dumux/freeflow/navierstokes/mass/1p/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Volume variables for the single-phase Navier-Stokes mass balance fused with the
 *        two-equation, high-Reynolds-number wall-function k-epsilon turbulence closure for the
 *        turbulent kinetic energy k and the dissipation rate epsilon (see turbulenceequations.md
 *        \S8 for the physics, and whatisimplemented.md/proposedimplementation.md for why k/epsilon
 *        live here rather than on a separate coupled sub-domain, following the same strategy
 *        already validated for the one-equation and k-omega models).
 *
 * Closure constants/functions ported near-verbatim from the deleted
 * releases/3.10:dumux/freeflow/rans/twoeq/kepsilon/volumevariables.hh.
 * stressTensorScalarProduct() is read from the problem (forwarded from the momentum domain's
 * Dumux::RANSMomentumProblem); inNearWallRegion()/isMatchingPoint()/yPlusNominal()/
 * uStarNominal() are read from the mass problem's own (lagged) per-time-step wall-function
 * bookkeeping (Dumux::KEpsilonMassProblem), which itself needs momentum-domain data (wall
 * distance, velocity gradients) forwarded the same way.
 */
template<class Traits>
class KEpsilonMassVolumeVariables
: public NavierStokesMassOnePVolumeVariables<Traits>
{
    using ParentType = NavierStokesMassOnePVolumeVariables<Traits>;
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    using Indices = typename Traits::ModelTraits::Indices;

    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        const auto eIdx = scv.elementIndex();
        stressTensorScalarProduct_ = problem.stressTensorScalarProduct(eIdx);
        inNearWallRegion_ = problem.inNearWallRegion(eIdx);
        isMatchingPoint_ = problem.isMatchingPoint(eIdx);
        yPlusNominal_ = problem.yPlusNominal(eIdx);
        uStarNominal_ = problem.uStarNominal(eIdx);

        if (inNearWallRegion_)
            dynamicEddyViscosity_ = problem.zeroEqDynamicEddyViscosity(eIdx);
        else
            dynamicEddyViscosity_ = cMu()*turbulentKineticEnergy()*turbulentKineticEnergy()/dissipation()*this->density();
    }

    //! The turbulent kinetic energy primary variable.
    Scalar turbulentKineticEnergy() const
    { return this->priVar(Indices::turbulentKineticEnergyIdx); }

    //! The dissipation rate primary variable.
    Scalar dissipation() const
    { return this->priVar(Indices::dissipationIdx); }

    Scalar stressTensorScalarProduct() const
    { return stressTensorScalarProduct_; }

    bool inNearWallRegion() const
    { return inNearWallRegion_; }

    bool isMatchingPoint() const
    { return isMatchingPoint_; }

    Scalar yPlusNominal() const
    { return yPlusNominal_; }

    Scalar uStarNominal() const
    { return uStarNominal_; }

    //! The dynamic eddy viscosity: the standard k-epsilon closure outside the near-wall region,
    //! the (lagged) blended zero-equation value inside it - the "two-layer model".
    Scalar dynamicEddyViscosity() const
    { return dynamicEddyViscosity_; }

    Scalar cMu() const { return 0.09; }
    Scalar sigmaK() const { return 1.0; }
    Scalar sigmaEpsilon() const { return 1.3; }
    Scalar cOneEpsilon() const { return 1.44; }
    Scalar cTwoEpsilon() const { return 1.92; }

protected:
    Scalar stressTensorScalarProduct_ = 0.0;
    Scalar yPlusNominal_ = 0.0;
    Scalar uStarNominal_ = 0.0;
    Scalar dynamicEddyViscosity_ = 0.0;
    bool inNearWallRegion_ = false;
    bool isMatchingPoint_ = false;
};

} // end namespace Dumux

#endif
