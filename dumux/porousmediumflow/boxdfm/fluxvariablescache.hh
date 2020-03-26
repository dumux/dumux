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
 * \ingroup BoxDFMModel
 * \brief Cache class for the flux variables to be used
 *        in conjunction with the box discrete fracture scheme.
 */

#ifndef DUMUX_POROUSMEDIUM_BOXDFM_FLUXVARIABLESCACHE_HH
#define DUMUX_POROUSMEDIUM_BOXDFM_FLUXVARIABLESCACHE_HH

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux {

/*!
 * \ingroup BoxDFMModel
 * \brief We only store discretization-related quantities for the box method.
 *        However, we cannot reuse the cache of the standard box method as we have
 *        to take into account the scvs that lie on fracture facets.
 */
template<class TypeTag>
class BoxDfmFluxVariablesCache
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridIndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<GridIndexType>;
    using TransmissibilityVector = std::vector<GridIndexType>;

    using CoordScalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
    using JacobianInverseTransposed = typename Element::Geometry::JacobianInverseTransposed;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:

    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace& scvf)
    {
        const auto geometry = element.geometry();
        const auto& localBasis = fvGeometry.feLocalBasis();

        // evaluate shape functions and gradients at the integration point
        std::vector<ShapeValue> shapeVals;
        const auto ipLocal = geometry.local(scvf.ipGlobal());
        jacInvT_ = geometry.jacobianInverseTransposed(ipLocal);
        localBasis.evaluateJacobian(ipLocal, shapeJacobian_);
        localBasis.evaluateFunction(ipLocal, shapeVals);

        // set the shape values
        shapeValues_.resize(fvGeometry.numScv(), 0.0);
        if (!scvf.isOnFracture())
            std::copy(shapeVals.begin(), shapeVals.end(), shapeValues_.begin());
        else
        {
            const auto thisFacetIdx = scvf.facetIndexInElement();
            for (const auto& scv: scvs(fvGeometry))
                if (scv.isOnFracture() && scv.facetIndexInElement() == thisFacetIdx)
                    shapeValues_[scv.indexInElement()] = shapeVals[scv.localDofIndex()];
        }

        // set the shape value gradients
        gradN_.resize(fvGeometry.numScv(), GlobalPosition(0.0));
        if (!scvf.isOnFracture())
        {
            for (const auto& scv: scvs(fvGeometry))
                if (!scv.isOnFracture())
                    jacInvT_.mv(shapeJacobian_[scv.localDofIndex()][0], gradN_[scv.indexInElement()]);
        }
        else
        {
            const auto thisFacetIdx = scvf.facetIndexInElement();

            // first, find all local dofs on this facet
            std::vector<unsigned int> facetLocalDofs;
            for (const auto& scv : scvs(fvGeometry))
                if (scv.isOnFracture() && scv.facetIndexInElement() == thisFacetIdx)
                    facetLocalDofs.push_back(scv.localDofIndex());

            for (const auto& scv: scvs(fvGeometry))
            {
                // now, create entries for all fracture scvs on this same facet ...
                if (scv.isOnFracture() && scv.facetIndexInElement() == thisFacetIdx)
                    jacInvT_.mv(shapeJacobian_[scv.localDofIndex()][0], gradN_[scv.indexInElement()]);

                // ... and those non-fracture scvs that are not on this facet
                else if (!scv.isOnFracture()
                         && std::find( facetLocalDofs.begin(),
                                       facetLocalDofs.end(),
                                       scv.localDofIndex() ) == facetLocalDofs.end())
                {
                    jacInvT_.mv(shapeJacobian_[scv.localDofIndex()][0], gradN_[scv.indexInElement()]);
                }
            }
        }
    }

    //! Returns the Jacobian of the shape functions at the integration point.
    const std::vector<ShapeJacobian>& shapeJacobian() const { return shapeJacobian_; }
    //! Returns the shape values for all scvs at the integration point.
    const std::vector<ShapeValue>& shapeValues() const { return shapeValues_; }
    //! Returns the shape value gradients for all scvs at the integration point.
    const JacobianInverseTransposed& jacInvT() const { return jacInvT_; }
    //! Returns the shape value gradients corresponding to an scv.
    const GlobalPosition& gradN(unsigned int scvIdxInElement) const { return gradN_[scvIdxInElement]; }

private:
    std::vector<GlobalPosition> gradN_;
    std::vector<ShapeJacobian> shapeJacobian_;
    std::vector<ShapeValue> shapeValues_;
    JacobianInverseTransposed jacInvT_;
};

} // end namespace Dumux

#endif
