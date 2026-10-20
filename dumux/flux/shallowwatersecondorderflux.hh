// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \copydoc Dumux::ShallowWaterSecondOrderFlux
 */
#ifndef DUMUX_FLUX_SHALLOW_WATER_SECOND_ORDER_FLUX_HH
#define DUMUX_FLUX_SHALLOW_WATER_SECOND_ORDER_FLUX_HH

#include <type_traits>
#include <utility>

#include <dune/common/exceptions.hh>

#include <dumux/discretization/cellcentered/tworingconnectivitymap.hh>
#include <dumux/flux/fluxvariablescaching.hh>
#include <dumux/flux/shallowwater/reconstruction.hh>
#include <dumux/flux/shallowwater/riemannproblem.hh>

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief Second-order shallow water advective flux with limited linear reconstruction.
 *
 * The states of both cells are reconstructed at the face (ShallowWater::CellReconstruction) and passed to
 * the same Riemann problem with hydrostatic reconstruction as in ShallowWaterFlux (ShallowWater::riemannProblem).
 * With a reconstructed bed, the hydrostatic reconstruction needs a centered bed slope source in each cell
 * (Audusse et al. 2004, \cite Audusse2004). It is split into contributions of the faces of the cell,
 * \f$ S_f = g h_i (z_f - z_i) - \frac{g}{2} (h_f - h_i)^2 \f$ times the unit outer normal and the face area,
 * where \f$ h_f \f$ and \f$ z_f \f$ are the values of cell \f$ i \f$ reconstructed at the face, and added
 * to the momentum flux of cell \f$ i \f$ over the face (for boundary faces see boundaryBedSlopeSource).
 * Together with the hydrostatic pressure \f$ \frac{g}{2} h_f^2 \f$ of the face fluxes, the sum over the
 * faces is \f$ g h_i |V_i| \nabla \eta \f$ for linear fields and vanishes for a lake at rest. Like the
 * fluxes, the contributions depend on the face neighbors, which the assembly of the cell-centered schemes
 * accounts for only in fluxes, not in sources. Boundary fluxes should be computed with the reconstructed
 * inside state (insideFaceState).
 *
 * The flux over a face depends on both cells and on all their face neighbors, so the grid geometry has
 * to use the two-ring connectivity map (CCTpfaTwoRingGridGeometryTraits), and the grid geometry and the
 * grid volume variables have to be cached. Distributed grids are not supported.
 *
 * \note For time steps at which the storage term does not dominate the Jacobian (wave Courant numbers much
 *       larger than one), the Jacobian is not diagonally dominant, and the iterative solvers with ILU or AMG
 *       preconditioners may not converge; a direct solver does.
 * \note With the limiter, the Newton method for the stationary equations may not converge; implicit Euler
 *       steps of growing size reach the steady state.
 */
template<class NumEqVector>
class ShallowWaterSecondOrderFlux
{
    using Scalar = typename NumEqVector::value_type;

public:
    using Cache = FluxVariablesCaching::EmptyAdvectionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller;

    /*!
     * \brief The advective flux over an inner face
     */
    template<class Problem, class FVElementGeometry, class ElementVolumeVariables>
    static NumEqVector flux(const Problem& problem,
                            const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {
        checkSetup_(fvGeometry);
        using Reconstruction = ShallowWater::CellReconstruction<Scalar, typename FVElementGeometry::SubControlVolumeFace::GlobalPosition>;
        const auto& params = parameters_(problem);

        const auto insideIdx = scvf.insideScvIdx();
        const auto outsideIdx = scvf.outsideScvIdx();
        const auto& nxy = scvf.unitOuterNormal();

        const auto gravity = problem.spatialParams().gravity(scvf.center());
        const Reconstruction inside(params, fvGeometry, elemVolVars, insideIdx, gravity);
        const Reconstruction outside(params, fvGeometry, elemVolVars, outsideIdx, gravity);
        const auto left = inside.faceState(scvf.center(), nxy, elemVolVars[outsideIdx], fvGeometry.scv(outsideIdx).center());
        const auto right = outside.faceState(scvf.center(), nxy, elemVolVars[insideIdx], fvGeometry.scv(insideIdx).center());

        const auto riemannFlux = ShallowWater::riemannProblem(left.waterDepth, right.waterDepth,
                                                              left.velocityX, right.velocityX,
                                                              left.velocityY, right.velocityY,
                                                              left.bedSurface, right.bedSurface,
                                                              gravity, nxy);

        NumEqVector localFlux(0.0);
        localFlux[0] = riemannFlux[0]*scvf.area();
        localFlux[1] = riemannFlux[1]*scvf.area();
        localFlux[2] = riemannFlux[2]*scvf.area();
        localFlux += bedSlopeSource_(inside, left, gravity, scvf);
        return localFlux;
    }

    /*!
     * \brief The state of the inside cell reconstructed at a boundary face
     */
    template<class Problem, class FVElementGeometry, class ElementVolumeVariables>
    static ShallowWater::FaceState<Scalar> insideFaceState(const Problem& problem,
                                                           const FVElementGeometry& fvGeometry,
                                                           const ElementVolumeVariables& elemVolVars,
                                                           const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {
        checkSetup_(fvGeometry);
        using Reconstruction = ShallowWater::CellReconstruction<Scalar, typename FVElementGeometry::SubControlVolumeFace::GlobalPosition>;
        const auto gravity = problem.spatialParams().gravity(scvf.center());
        const Reconstruction inside(parameters_(problem), fvGeometry, elemVolVars, scvf.insideScvIdx(), gravity);
        return inside.faceState(scvf.center(), scvf.unitOuterNormal());
    }

    /*!
     * \brief The contribution of a boundary face to the centered bed slope source of the inside cell,
     *        to be added to the boundary flux
     */
    template<class Problem, class FVElementGeometry, class ElementVolumeVariables>
    static NumEqVector boundaryBedSlopeSource(const Problem& problem,
                                              const FVElementGeometry& fvGeometry,
                                              const ElementVolumeVariables& elemVolVars,
                                              const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {
        checkSetup_(fvGeometry);
        using Reconstruction = ShallowWater::CellReconstruction<Scalar, typename FVElementGeometry::SubControlVolumeFace::GlobalPosition>;
        const auto gravity = problem.spatialParams().gravity(scvf.center());
        const Reconstruction inside(parameters_(problem), fvGeometry, elemVolVars, scvf.insideScvIdx(), gravity);
        return bedSlopeSource_(inside, inside.faceState(scvf.center(), scvf.unitOuterNormal()), gravity, scvf);
    }

private:
    template<class Reconstruction, class SubControlVolumeFace>
    static NumEqVector bedSlopeSource_(const Reconstruction& cell,
                                       const ShallowWater::FaceState<Scalar>& state,
                                       const Scalar gravity,
                                       const SubControlVolumeFace& scvf)
    {
        const Scalar deltaDepth = state.waterDepth - cell.waterDepth();
        const Scalar faceSource = gravity*cell.waterDepth()*(state.bedSurface - cell.bedSurface())
                                  - 0.5*gravity*deltaDepth*deltaDepth;
        NumEqVector source(0.0);
        source[1] = faceSource*scvf.unitOuterNormal()[0]*scvf.area();
        source[2] = faceSource*scvf.unitOuterNormal()[1]*scvf.area();
        return source;
    }

    template<class Problem>
    static const ShallowWater::ReconstructionParameters<Scalar>& parameters_(const Problem& problem)
    {
        static const auto params = ShallowWater::readReconstructionParameters<Scalar>(problem.paramGroup());
        return params;
    }

    template<class FVElementGeometry>
    static void checkSetup_(const FVElementGeometry& fvGeometry)
    {
        using GridGeometry = typename FVElementGeometry::GridGeometry;
        using ConnectivityMap = std::decay_t<decltype(std::declval<const GridGeometry&>().connectivityMap())>;
        static_assert(std::is_same_v<ConnectivityMap, CCTwoRingConnectivityMap<GridGeometry>>,
                      "The second-order shallow water flux requires a grid geometry with the two-ring connectivity map "
                      "(CCTpfaTwoRingGridGeometryTraits)");

        static const bool distributed = fvGeometry.gridGeometry().gridView().comm().size() > 1;
        if (distributed)
            DUNE_THROW(Dune::NotImplemented, "The second-order shallow water flux on distributed grids");
    }
};

} // end namespace Dumux

#endif
