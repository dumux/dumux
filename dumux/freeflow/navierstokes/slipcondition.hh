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

template<class GridGeometry, class SlipCondition>
class SlipVelocityHelper;

/*!
 * \ingroup NavierStokesModel
 * \brief Navier Stokes slip velocity helper for Beavers-Joseph condition
 */
template<class GridGeometry>
class SlipVelocityHelper<GridGeometry, SlipConditions::BJ>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    static constexpr int dimWorld = GridView::dimensionworld;
    using Scalar = typename GridView::ctype;
    using Vector = Dune::FieldVector<Scalar, dimWorld>;

public:
    template<class Problem, class ElementVolumeVariables>
    static Vector velocity(const Problem& problem,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvf,
                           const ElementVolumeVariables& elemVolVars,
                           Scalar tangentialVelocityGradient)
   {
        if constexpr (GridGeometry::discMethod == DiscretizationMethods::fcstaggered)
        {
            assert(scvf.isLateral());
            assert(scvf.boundary());

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
                tangentialVelocityGradient = 0.0;

            const auto& porousMediumVelocity = problem.porousMediumVelocity(fvGeometry, scvf);
            const Scalar scalarSlipVelocity = (tangentialVelocityGradient*distanceNormalToBoundary
                + porousMediumVelocity * tangent * betaBJ * distanceNormalToBoundary
                + elemVolVars[scv].velocity()) / (betaBJ*distanceNormalToBoundary + 1.0);

            return scalarSlipVelocity*tangent;
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Slip velocity currently only implemented for fcstaggered");
    }
};

/*!
 * \ingroup NavierStokesModel
 * \brief Navier Stokes slip velocity helper for Beavers-Joseph-Saffman condition
 */
template<class GridGeometry>
class SlipVelocityHelper<GridGeometry, SlipConditions::BJS>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    static constexpr int dimWorld = GridView::dimensionworld;
    using Scalar = typename GridView::ctype;
    using Vector = Dune::FieldVector<Scalar, dimWorld>;

public:
    template<class Problem, class ElementVolumeVariables>
    static Vector velocity(const Problem& problem,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvf,
                           const ElementVolumeVariables& elemVolVars,
                           Scalar tangentialVelocityGradient)
   {
        if constexpr (GridGeometry::discMethod == DiscretizationMethods::fcstaggered)
        {
            assert(scvf.isLateral());
            assert(scvf.boundary());

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
                tangentialVelocityGradient = 0.0;
            const Scalar scalarSlipVelocity = (tangentialVelocityGradient*distanceNormalToBoundary
                + elemVolVars[scv].velocity()) / (betaBJ*distanceNormalToBoundary + 1.0);

            return scalarSlipVelocity*tangent;
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Slip velocity currently only implemented for fcstaggered");
    }
};


} // end namespace Dumux

#endif
