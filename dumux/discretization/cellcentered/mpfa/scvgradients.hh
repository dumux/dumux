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
 * \ingroup CCMpfaDiscretization
 * \brief Class roviding functionality for the reconstruction of the
 *        gradients in the sub-control volumes involved in mpfa schemes.
 */
#ifndef DUMUX_CC_MPFA_SCV_GRADIENTS_HH
#define DUMUX_CC_MPFA_SCV_GRADIENTS_HH

#include <vector>
#include <utility>
#include <type_traits>

#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/discretization/cellcentered/mpfa/localassemblerhelper.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class roviding functionality for the reconstruction of the
 *        gradients in the sub-control volumes involved in mpfa schemes.
 */
class CCMpfaScvGradients
{
    //! Return type of the gradient computation function (pair of scv centers and corresponding gradients)
    template<class GridGeometry, class Scalar>
    using ResultPair = std::pair< std::vector<typename GridGeometry::SubControlVolume::GlobalPosition>,
                                  std::vector<Dune::FieldVector<Scalar, GridGeometry::GridView::dimension>> >;

public:

    /*!
     * \brief Computes the phase velocities in the scvs of the grid.
     *
     * \param gridGeometry The finite volume grid geometry
     * \param gridVariables The variables living on the grid
     * \param x The vector containing the solution
     * \param phaseIdx The index of the fluid phase to be considered
     */
    template<class GridGeometry, class GridVariables, class SolutionVector>
    static ResultPair<GridGeometry, typename GridVariables::Scalar>
    computeVelocities(const GridGeometry& gridGeometry,
                      const GridVariables& gridVariables,
                      const SolutionVector& x,
                      unsigned int phaseIdx)
    {
        auto gradToVelocity = [phaseIdx] (const auto& grad, const auto& volVars)
        {
            auto vel = mv(volVars.permeability(), grad);
            vel *= volVars.mobility(phaseIdx);
            return vel;
        };

        return computePressureGradients(gridGeometry, gridVariables, x, phaseIdx, gradToVelocity);
    }

    /*!
     * \brief Computes the pressure gradients in the scvs of the grid.
     *
     * \param gridGeometry The finite volume grid geometry
     * \param gridVariables The variables living on the grid
     * \param x The vector containing the solution
     * \param phaseIdx The index of the fluid phase to be considered
     */
    template<class GridGeometry, class GridVariables, class SolutionVector>
    static ResultPair<GridGeometry, typename GridVariables::Scalar>
    computePressureGradients(const GridGeometry& gridGeometry,
                             const GridVariables& gridVariables,
                             const SolutionVector& x,
                             unsigned int phaseIdx)
    {
        auto f = [] (const auto& grad, const auto& volVars) { return grad; };
        return computePressureGradients(gridGeometry, gridVariables, x, phaseIdx, f);
    }

    /*!
     * \brief Computes the pressure gradients in the scvs of the grid.
     *
     * \param gridGeometry The finite volume grid geometry
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
    template<class GridGeometry, class GridVariables, class SolutionVector, class F>
    static ResultPair<GridGeometry, typename GridVariables::Scalar>
    computePressureGradients(const GridGeometry& gridGeometry,
                             const GridVariables& gridVariables,
                             const SolutionVector& x,
                             unsigned int phaseIdx,
                             F& f)
    {
        using ElemVolVars = typename GridVariables::GridVolumeVariables::LocalView;
        using FVElementGeometry = typename GridGeometry::LocalView;

        auto handleFunction = [&] (auto& result,
                                   const auto& handle,
                                   const auto& iv,
                                   const FVElementGeometry& fvGeometry,
                                   const ElemVolVars& elemVolVars)
        {
            handle.advectionHandle().setPhaseIndex(phaseIdx);
            const auto grads = InteractionVolumeAssemblerHelper::assembleScvGradients(handle.advectionHandle(), iv);
            for (unsigned int i = 0; i < iv.numScvs(); ++i)
            {
                const auto& volVars = elemVolVars[iv.localScv(i).gridScvIndex()];
                const auto scvGeometry = iv.getScvGeometry(i, fvGeometry);
                result.first.push_back(scvGeometry.center());
                result.second.push_back( f(grads[i], volVars) );
            }
        };

        return computeGradients_(gridGeometry, gridVariables, x, handleFunction);
    }

private:
    /*!
     * \brief Computes the gradients executing the provided function
     *        on an interaction volume / data handle pair.
     */
    template<class GridGeometry, class GridVariables, class SolutionVector, class HandleFunction>
    static ResultPair<GridGeometry, typename GridVariables::Scalar>
    computeGradients_(const GridGeometry& gridGeometry,
                      const GridVariables& gridVariables,
                      const SolutionVector& x,
                      const HandleFunction& handleFunction)
    {
        using GridView = typename GridGeometry::GridView;
        static constexpr int dim = GridView::dimension;

        // first, find out how many scvs live on this grid
        std::size_t numScvs = 0;
        const auto& gridView = gridGeometry.gridView();
        for (const auto& element : elements(gridView))
            numScvs += element.subEntities(dim);

        ResultPair<GridGeometry, typename GridVariables::Scalar> result;
        result.first.reserve(numScvs);
        result.second.reserve(numScvs);
        std::vector<bool> vertexHandled(gridView.size(dim), false);

        for (const auto& element : elements(gridView))
        {
            bool allFinished = true;
            for (int i = 0; i < element.subEntities(dim); ++i)
                if (!vertexHandled[gridGeometry.vertexMapper().subIndex(element, i, dim)])
                    allFinished = false;

            // bind views only if there is unfinished buisness
            if (allFinished)
                continue;

            // compute gradients in all scvs of all interaction volumes in this element
            auto fvGeometry = localView(gridGeometry);
            auto elemVolVars = localView(gridVariables.curGridVolVars());
            auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, x);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (vertexHandled[scvf.vertexIndex()])
                    continue;

                if (gridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
                {
                    const auto& iv = elemFluxVarsCache.secondaryInteractionVolume(scvf);
                    const auto& handle = elemFluxVarsCache.secondaryDataHandle(scvf);
                    handleFunction(result, handle, iv, fvGeometry, elemVolVars);
                }
                else
                {
                    const auto& iv = elemFluxVarsCache.primaryInteractionVolume(scvf);
                    const auto& handle = elemFluxVarsCache.primaryDataHandle(scvf);
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
