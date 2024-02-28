// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief Navier Stokes slip condition
 */
#ifndef DUMUX_NAVIERSTOKES_SLIPCONDITION_HH
#define DUMUX_NAVIERSTOKES_SLIPCONDITION_HH

#include <dumux/discretization/method.hh>

namespace Dumux {

namespace SlipConditions {

struct BJ : public Utility::Tag<BJ> {
    static std::string name() { return "Beavers-Joseph"; }
};

struct BJS : public Utility::Tag<BJS> {
    static std::string name() { return "Beavers-Joseph-Saffman"; }
};

inline constexpr BJ bj{};
inline constexpr BJS bjs{};

} // end namespace SlipModel

template<class DiscretizationMethod, class SlipCondition>
class SlipVelocityHelper;

/*!
 * \ingroup NavierStokesModel
 * \brief Navier Stokes slip velocity helper for fcstaggered discretization
 */
template< class SlipCondition>
struct SlipVelocityHelper<DiscretizationMethods::FCStaggered, SlipCondition>
{
    static constexpr SlipCondition slipCondition{};

    /*!
     * \brief Returns the slip velocity at a porous boundary based on the Beavers-Joseph(-Saffman) condition.
     * \note This only returns a vector filled with one component of the slip velocity (corresponding to the dof axis of the scv the svf belongs to)
     * \param problem The problem
     * \param fvGeometry The finite-volume geometry
     * \param scvf The sub control volume face
     * \param elemVolVars The volume variables for the element
     * \param tangentialVelocityDeriv Pre-calculated velocity derivative
     */
    template<class Problem, class FVElementGeometry, class ElementVolumeVariables, class Scalar>
    static auto velocity(const Problem& problem,
                           const FVElementGeometry& fvGeometry,
                           const typename FVElementGeometry::SubControlVolumeFace& scvf,
                           const ElementVolumeVariables& elemVolVars,
                           Scalar tangentialVelocityDeriv)
   {
        assert(scvf.isLateral());
        assert(scvf.boundary());

        static constexpr int dimWorld = FVElementGeometry::GridGeometry::GridView::dimensionworld;
        using Vector = Dune::FieldVector<Scalar, dimWorld>;

        Vector porousMediumVelocity(0.0);

        if constexpr (slipCondition == SlipConditions::bj)
            porousMediumVelocity = problem.porousMediumVelocity(fvGeometry, scvf);
        else if (!(slipCondition == SlipConditions::bjs))
            DUNE_THROW(Dune::NotImplemented, "Fcstaggered currently only implements BJ or BJS slip conditions");

        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

        // create a unit normal vector oriented in positive coordinate direction
        Vector tangent(0.0);
        tangent[scv.dofAxis()] = 1.0;

        // du/dy + dv/dx = beta * (u_boundary-uPM)
        // beta = alpha/sqrt(K)
        const Scalar betaBJ = problem.betaBJ(fvGeometry, scvf, tangent);
        const Scalar distanceNormalToBoundary = (scvf.ipGlobal() - scv.dofPosition()).two_norm();

        static const bool onlyNormalGradient = getParamFromGroup<bool>(problem.paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph", false);
        if (onlyNormalGradient)
            tangentialVelocityDeriv = 0.0;

        const Scalar scalarSlipVelocity = (tangentialVelocityDeriv*distanceNormalToBoundary
            + porousMediumVelocity * tangent * betaBJ * distanceNormalToBoundary
            + elemVolVars[scv].velocity()) / (betaBJ*distanceNormalToBoundary + 1.0);

        return scalarSlipVelocity*tangent;
    }
};

} // end namespace Dumux

#endif
