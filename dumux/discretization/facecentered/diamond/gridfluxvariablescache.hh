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
 * \brief DOC ME
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_GRID_FLUXVARSCACHE_HH

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/facecentered/diamond/elementfluxvariablescache.hh>

namespace Dumux {

/*!
 * \brief DOC ME
 */
template<class P, class GG, class FVC, class FVCF>
struct FaceCenteredDiamondDefaultGridFVCTraits
{
    using Problem = P;
    using GridGeometry = GG;
    using FluxVariablesCache = FVC;
    using FluxVariablesCacheFiller = FVCF;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = FaceCenteredDiamondElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Flux variable caches on a gridview
 * \note The class is specialized for a version with and without grid caching
 */
template<class Problem,
         class GridGeometry,
         class FluxVariablesCache,
         class FluxVariablesCacheFiller,
         bool cachingEnabled = false,
         class Traits = FaceCenteredDiamondDefaultGridFVCTraits<Problem, GridGeometry, FluxVariablesCache, FluxVariablesCacheFiller>>
class FaceCenteredDiamondGridFluxVariablesCache;

//! Element-wise flux cache (variables that are the same for all scvfs in an element)
template<class Tensor>
struct FluxVariablesElementCache
{
    Tensor avgSkewGradV;
};

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Flux variable caches on a gridview with grid caching enabled
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 */
template<class P, class GG, class FVC, class FVCF, class Traits>
class FaceCenteredDiamondGridFluxVariablesCache<P, GG, FVC, FVCF, true, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = FaceCenteredDiamondGridFluxVariablesCache<P, GG, FVC, FVCF, true, Traits>;
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

    // When global caching is enabled, precompute transmissibilities and stencils for all the scv faces
    template<class GridVolumeVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false)
    {
        // only do the update if fluxes are solution dependent or if update is forced
        fluxVarsCache_.resize(gridGeometry.numScvf());
        fluxVarsElementCache_.resize(gridGeometry.gridView().size(0));
        fluxVarsElementCenterGradN_.resize(gridGeometry.gridView().size(0));
        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            // Prepare the geometries within the elements of the stencil
            fvGeometry.bind(element);

            auto elemVolVars = localView(gridVolVars);
            elemVolVars.bind(element, fvGeometry, sol);

            const auto eIdx = gridGeometry.elementMapper().index(element);

            // only update shape functions at element centers if update is forced
            if (forceUpdate)
                fluxVarsElementCenterGradN_[eIdx] = computeCenterShapeFunctionGradients_(fvGeometry);

            // compute element-wise tensor
            fluxVarsElementCache_[eIdx] = FluxVariablesElementCache<Tensor>{
                neighborhoodAverageSkewGradV_(fvGeometry, gridVolVars)
            };

            // only update shape functions for fluxes if update is forced
            if (forceUpdate)
                for (const auto& scvf : scvfs(fvGeometry))
                    fluxVarsCache_[scvf.index()].update(problem(), element, fvGeometry, elemVolVars, scvf);
        }
    }

    template<class ElementVolumeVariables>
    void updateElement(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars)
    {
        // update element-wise tensor
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        fluxVarsElementCache_[eIdx] = FluxVariablesElementCache<Tensor>{
            neighborhoodAverageSkewGradV_(fvGeometry, elemVolVars.gridVolVars())
        };
    }

    const Problem& problem() const
    { return *problemPtr_; }

    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.index()]; }

    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.index()]; }

    const FluxVariablesElementCache<Tensor>& operator [](std::size_t eIdx) const
    { return fluxVarsElementCache_[eIdx]; }

    FluxVariablesElementCache<Tensor>& operator [](std::size_t eIdx)
    { return fluxVarsElementCache_[eIdx]; }

private:
    std::vector<GlobalPosition> computeCenterShapeFunctionGradients_(const FVElementGeometry& fvGeometry) const
    {
        const auto& element = fvGeometry.element();
        const auto& geometry = element.geometry();
        const auto ipLocal = geometry.local(geometry.center());

        const auto& localBasis = fvGeometry.feLocalBasis();

        // TODO cache these values
        const auto jacInvT = geometry.jacobianInverseTransposed(ipLocal);
        std::vector<ShapeJacobian> shapeJacobian;
        localBasis.evaluateJacobian(ipLocal, shapeJacobian);

        // compute the gradN at for every scv/dof
        std::vector<GlobalPosition> gradN(fvGeometry.numScv(), GlobalPosition(0));
        for (const auto& scv : scvs(fvGeometry))
            jacInvT.mv(shapeJacobian[scv.localDofIndex()][0], gradN[scv.indexInElement()]);

        return gradN;
    }

    template<class GridVolumeVariables>
    Tensor elementAverageSkewGradV_(const FVElementGeometry& fvGeometry, const GridVolumeVariables& gridVolVars) const
    {
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(fvGeometry.element());
        const auto& gradN = fluxVarsElementCenterGradN_[eIdx];

        // integrate tensor with midpoint rule
        Tensor skewGradV(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& thisGradN = gradN[scv.indexInElement()];
            const auto& volVars = gridVolVars.volVars(scv);
            for (int dir = 0; dir < dim; ++dir)
                skewGradV[dir].axpy(0.5*volVars.velocity(dir), thisGradN);
        }

        // we need the skew-symmetric part here, i.e. 0.5*(gradV - gradV^T)
        skewGradV.axpy(-0.5, getTransposed(skewGradV));

        return skewGradV;
    }

    template<class GridVolumeVariables>
    Tensor neighborhoodAverageSkewGradV_(const FVElementGeometry& fvGeometry, const GridVolumeVariables& gridVolVars) const
    {
        const auto& element = fvGeometry.element();
        Scalar totalVolume = element.geometry().volume();
        Tensor skewGradV(0.0);
        skewGradV += elementAverageSkewGradV_(fvGeometry, gridVolVars) * totalVolume;
        for (const auto& intersection : intersections(fvGeometry.gridGeometry().gridView(), element))
        {
            if (intersection.neighbor())
            {
                const auto neighbor = intersection.outside();
                const auto fvGeometryNeighbor = localView(fvGeometry.gridGeometry()).bind(neighbor);
                const auto vol = neighbor.geometry().volume();
                skewGradV += elementAverageSkewGradV_(fvGeometry, gridVolVars) * vol;
                totalVolume += vol;
            }
        }

        skewGradV /= totalVolume;

        return skewGradV;
    }

    const Problem* problemPtr_;
    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<FluxVariablesElementCache<Tensor>> fluxVarsElementCache_;
    std::vector<std::vector<GlobalPosition>> fluxVarsElementCenterGradN_;
};

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Flux variable caches on a gridview with grid caching disabled
 */
template<class P, class GG, class FVC, class FVCF, class Traits>
class FaceCenteredDiamondGridFluxVariablesCache<P, GG, FVC, FVCF, false, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = FaceCenteredDiamondGridFluxVariablesCache<P, GG, FVC, FVCF, false, Traits>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    FaceCenteredDiamondGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false) {}

    template<class FVElementGeometry, class ElementVolumeVariables>
    void updateElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars) {}

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif