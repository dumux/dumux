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
 * \ingroup CCMpfaDiscretization
 * \brief Class roviding functionality for the reconstruction of the
 *        gradients in the sub-control volumes involved in mpfa schemes.
 */
#ifndef DUMUX_CC_MPFA_SCV_GRADIENTS_HH
#define DUMUX_CC_MPFA_SCV_GRADIENTS_HH

#include <type_traits>

#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class roviding functionality for the reconstruction of the
 *        gradients in the sub-control volumes involved in mpfa schemes.
 */
class CCMpfaScvGradients
{
    //! Return type of the gradient computation function
    template<class FVGridGeometry, class Scalar>
    struct ScvGradients
    {
        using GlobalPosition = typename FVGridGeometry::SubControlVolume::GlobalPosition;
        using Gradient = Dune::FieldVector<Scalar, FVGridGeometry::GridView::dimension>;

        std::vector<GlobalPosition> scvCenters;
        std::vector<Gradient> scvGradients;
    };

public:

    /*!
     * \brief Computes the phase velocities in the scvs of the grid.
     *
     * \param fvGridGeometry The finite volume grid geometry
     * \param gridVariables The variables living on the grid
     * \param x The vector containing the solution
     * \param phaseIdx The index of the fluid phase to be considered
     */
    template<class FVGridGeometry, class GridVariables, class SolutionVector>
    static ScvGradients<FVGridGeometry, typename GridVariables::Scalar>
    computeVelocities(const FVGridGeometry& fvGridGeometry,
                      const GridVariables& gridVariables,
                      const SolutionVector& x,
                      unsigned int phaseIdx)
    {
        using ResultType = ScvGradients<FVGridGeometry, typename GridVariables::Scalar>;
        using Gradient = typename ResultType::Gradient;

        auto gradToVelocity = [phaseIdx] (const Gradient& grad, const auto& volVars)
        {
            auto vel = mv(volVars.permeability(), grad);
            vel *= volVars.mobility(phaseIdx);
            return vel;
        };

        return computePressureGradients(fvGridGeometry, gridVariables, x, phaseIdx, gradToVelocity);
    }

    /*!
     * \brief Computes the pressure gradients in the scvs of the grid.
     *
     * \param fvGridGeometry The finite volume grid geometry
     * \param gridVariables The variables living on the grid
     * \param x The vector containing the solution
     * \param phaseIdx The index of the fluid phase to be considered
     */
    template<class FVGridGeometry, class GridVariables, class SolutionVector>
    static ScvGradients<FVGridGeometry, typename GridVariables::Scalar>
    computePressureGradients(const FVGridGeometry& fvGridGeometry,
                             const GridVariables& gridVariables,
                             const SolutionVector& x,
                             unsigned int phaseIdx)
    {
        auto f = [] (const auto& grad, const auto& volVars) { return grad; };
        return computePressureGradients(fvGridGeometry, gridVariables, x, phaseIdx, f);
    }

    /*!
     * \brief Computes the pressure gradients in the scvs of the grid.
     *
     * \param fvGridGeometry The finite volume grid geometry
     * \param gridVariables The variables living on the grid
     * \param x The vector containing the solution
     * \param phaseIdx The index of the fluid phase to be considered
     * \param f a function which is applied to the gradients and which
     *          has the following signature:
     *
     *          Gradient f(const Gradient& gradient, const VolumeVariables& volVars)
     *
     *          It receives the scv-gradient and the corresponding volume
     *          variables and returns the modified gradient. This can be
     *          used e.g. to turn the pressure gradients into velocities.
     */
    template<class FVGridGeometry, class GridVariables, class SolutionVector, class F>
    static ScvGradients<FVGridGeometry, typename GridVariables::Scalar>
    computePressureGradients(const FVGridGeometry& fvGridGeometry,
                             const GridVariables& gridVariables,
                             const SolutionVector& x,
                             unsigned int phaseIdx,
                             F& f)
    {
        using ResultType = ScvGradients<FVGridGeometry, typename GridVariables::Scalar>;
        using Problem = std::decay_t<decltype(gridVariables.curGridVolVars().problem())>;
        using ElemVolVars = typename GridVariables::GridVolumeVariables::LocalView;
        using FVElementGeometry = typename FVGridGeometry::LocalView;

        auto handleFunction = [&] (ResultType& result,
                                   const auto& handle,
                                   const auto& iv,
                                   const FVElementGeometry& fvGeometry,
                                   const ElemVolVars& elemVolVars)
        {
            using IV = std::decay_t<decltype(iv)>;
            using LocalAssembler = typename IV::Traits::template LocalAssembler<Problem, FVElementGeometry, ElemVolVars>;

            handle.advectionHandle().setPhaseIndex(phaseIdx);
            const auto grads = LocalAssembler::assembleScvGradients(handle.advectionHandle(), iv);
            for (unsigned int i = 0; i < iv.numScvs(); ++i)
            {
                const auto& volVars = elemVolVars[iv.localScv(i).gridScvIndex()];
                const auto scvGeometry = iv.getScvGeometry(i, fvGeometry);
                result.scvCenters.push_back(scvGeometry.center());
                result.scvGradients.push_back( f(grads[i], volVars) );
            }
        };

        return computeGradients_(fvGridGeometry, gridVariables, x, handleFunction);
    }

private:
    /*!
     * \brief Computes the gradients executing the provided function
     *        on an interaction volume / data handle pair.
     */
    template<class FVGridGeometry, class GridVariables, class SolutionVector, class HandleFunction>
    static ScvGradients<FVGridGeometry, typename GridVariables::Scalar>
    computeGradients_(const FVGridGeometry& fvGridGeometry,
                      const GridVariables& gridVariables,
                      const SolutionVector& x,
                      const HandleFunction& handleFunction)
    {
        using GridView = typename FVGridGeometry::GridView;
        static constexpr int dim = GridView::dimension;

        // first, find out how many scvs live on this grid
        std::size_t numScvs = 0;
        const auto& gridView = fvGridGeometry.gridView();
        for (const auto& element : elements(gridView))
            numScvs += element.subEntities(dim);

        ScvGradients<FVGridGeometry, typename GridVariables::Scalar> result;
        result.scvCenters.reserve(numScvs);
        result.scvGradients.reserve(numScvs);
        std::vector<bool> vertexHandled(gridView.size(dim), false);

        for (const auto& element : elements(gridView))
        {
            bool allFinished = true;
            for (int i = 0; i < element.subEntities(dim); ++i)
                if (!vertexHandled[fvGridGeometry.vertexMapper().subIndex(element, i, dim)])
                    allFinished = false;

            // bind views only if there is unfinished buisness
            if (allFinished)
                continue;

            // compute gradients in all scvs of all interaction volumes in this element
            auto fvGeometry = localView(fvGridGeometry);
            auto elemVolVars = localView(gridVariables.curGridVolVars());
            auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, x);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (vertexHandled[scvf.vertexIndex()])
                    continue;

                const auto& fluxVarsCache = elemFluxVarsCache[scvf];
                if (fvGridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
                {
                    const auto& iv = elemFluxVarsCache.secondaryInteractionVolumes()[fluxVarsCache.ivIndexInContainer()];
                    const auto& handle = elemFluxVarsCache.secondaryDataHandles()[fluxVarsCache.ivIndexInContainer()];
                    handleFunction(result, handle, iv, fvGeometry, elemVolVars);
                }
                else
                {
                    const auto& iv = elemFluxVarsCache.primaryInteractionVolumes()[fluxVarsCache.ivIndexInContainer()];
                    const auto& handle = elemFluxVarsCache.primaryDataHandles()[fluxVarsCache.ivIndexInContainer()];
                    handleFunction(result, handle, iv, fvGeometry, elemVolVars);
                }

                vertexHandled[scvf.vertexIndex()] = true;
            }
        }

        return result;
    }
};

} // end namespace Dumux

#endif
