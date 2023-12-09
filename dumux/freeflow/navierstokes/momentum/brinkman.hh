// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief A helper function to add Brinkman term to the momentum balance
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_BRINKMAN_TERM_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_BRINKMAN_TERM_HH

#include <type_traits>
#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/math.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/freeflow/navierstokes/momentum/velocityreconstruction.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief A helper function to add Brinkman term to the momentum balance
 * \addtogroup NavierStokesModel
 * @{
 * The Navier-Stokes model can be extended to a Darcy-Brinkman model by adding
 * the term:
 * \f[
      + \epsilon_B \mu \mathbf{K}^{-1} \mathbf{v}
 * \f]
 * to the momentum balance. This can be achieved with the helper function \ref addBrinkmanTerm.
 * The function relies on the spatial parameters class being based on <code>BrinkmanSpatialParams</code>
 * or providing the `brinkmanEpsilon` and `inversePermeability` interfaces.
 * These interface functions provide the
 * weighting factor \f$ \epsilon_B \f$ and the permeability tensor \f$ \mathbf{K} \f$.
 * @}
 */
template<class NumEqVector, class Problem, class FVElementGeometry, class ElementVolumeVariables>
void addBrinkmanTerm(
    NumEqVector& source,
    const Problem& problem,
    const typename FVElementGeometry::Element& element,
    const FVElementGeometry& fvGeometry,
    const ElementVolumeVariables& elemVolVars,
    const typename FVElementGeometry::SubControlVolume& scv
){
    if (problem.spatialParams().brinkmanEpsilon(element, fvGeometry, scv) > 0.0)
    {
        const auto brinkmanEpsilon = problem.spatialParams().brinkmanEpsilon(element, fvGeometry, scv);
        const auto invK = problem.spatialParams().inversePermeability(element, fvGeometry, scv);
        const auto mu = problem.effectiveViscosity(element, fvGeometry, scv);

        using DM = typename FVElementGeometry::GridGeometry::DiscretizationMethod;
        if constexpr (DiscretizationMethods::isCVFE<DM>)
        {
            const auto velocity = elemVolVars[scv].velocity();
            source -= mu * brinkmanEpsilon * mv(invK, velocity);
        }
        else if constexpr (DM{} == DiscretizationMethods::fcstaggered)
        {
            if constexpr (Dune::IsNumber<std::decay_t<decltype(invK)>>::value)
            {
                const auto velocity = elemVolVars[scv].velocity();
                source -= mu * brinkmanEpsilon * invK * velocity;
            }
            else
            {
                // permeability is tensor-valued, use velocity vector reconstruction
                const auto getVelocitySCV = [&](const auto& scv){ return elemVolVars[scv].velocity(); };
                const auto velocity = StaggeredVelocityReconstruction::faceVelocity(scv, fvGeometry, getVelocitySCV);
                source -= mu * brinkmanEpsilon * mv(invK, velocity)[scv.dofAxis()];
            }
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Brinkman term only implemented for CVFE and Staggered");
    }
}

} // end namespace Dumux

#endif
