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
 * \ingroup DiamondDiscretization
 * \copydoc Dumux::FaceCenteredDiamondGridFluxVariablesCache
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_GRID_FLUXVARSCACHE_HH

#include <dumux/parallel/parallel_for.hh>

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cvfe/elementfluxvariablescache.hh>

namespace Dumux {

/*!
 * \ingroup DiamondDiscretization
 * \brief Default traits for the diamond discretization
 */
template<class P, class GG, class FVC, class FVCF>
struct FaceCenteredDiamondDefaultGridFVCTraits
{
    using Problem = P;
    using GridGeometry = GG;
    using FluxVariablesCache = FVC;
    using FluxVariablesCacheFiller = FVCF;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = CVFEElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup DiamondDiscretization
 * \brief Flux variable caches on a gridview
 */
template<class P, class GG, class FVC, class FVCF, bool enableCaching,
         class Traits = FaceCenteredDiamondDefaultGridFVCTraits<P, GG, FVC, FVCF>>
class FaceCenteredDiamondGridFluxVariablesCache
{
    using Problem = typename Traits::Problem;
    using ThisType = FaceCenteredDiamondGridFluxVariablesCache<P, GG, FVC, FVCF, enableCaching, Traits>;
    using GridGeometry = typename Traits::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Element = typename GridGeometry::GridView::Grid::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr int dim = GridGeometry::GridView::dimension;
    using FeLocalBasis = typename GridGeometry::FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;
    using ShapeValue = typename Dune::FieldVector<double, 1>;
    using Tensor = Dune::FieldMatrix<double, dim>;
    using Scalar = double;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! export the sub control volume face type
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    FaceCenteredDiamondGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridVolumeVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false)
    {
        // only do the update if fluxes are solution dependent or if update is forced
        if (forceUpdate)
        {
            fluxVarsCache_.resize(gridGeometry.gridView().size(0));
            Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
            {
                // Prepare the geometries within the elements of the stencil
                const auto element = gridGeometry.element(eIdx);
                const auto fvGeometry = localView(gridGeometry).bind(element);
                const auto elemVolVars = localView(gridVolVars).bind(element, fvGeometry, sol);

                // only update shape functions for fluxes if update is forced
                fluxVarsCache_[eIdx].resize(fvGeometry.numScvf());
                for (const auto& scvf : scvfs(fvGeometry))
                    cache(eIdx, scvf.index()).update(problem, element, fvGeometry, elemVolVars, scvf);
            });
        }
    }

    const Problem& problem() const
    { return *problemPtr_; }

    // access operator
    const FluxVariablesCache& cache(std::size_t eIdx, std::size_t scvfIdx) const
    { return fluxVarsCache_[eIdx][scvfIdx]; }

    // access operator
    FluxVariablesCache& cache(std::size_t eIdx, std::size_t scvfIdx)
    { return fluxVarsCache_[eIdx][scvfIdx]; }

private:
    // currently bound element
    const Problem* problemPtr_;
    std::vector<std::vector<FluxVariablesCache>> fluxVarsCache_;
};

} // end namespace Dumux

#endif
