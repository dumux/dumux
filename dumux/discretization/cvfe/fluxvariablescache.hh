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
 * \ingroup CvfeDiscretization
 * \brief Flux variables cache class for the cvfe scheme
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_FLUXVARIABLES_CACHE_HH
#define DUMUX_DISCRETIZATION_CVFE_FLUXVARIABLES_CACHE_HH

#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup CvfeDiscretization
 * \brief Flux variables cache class for the cvfe scheme.
 *        For the cvfe scheme, this class does not contain any physics-/process-dependent
 *        data. It solely stores disretization-/grid-related data.
 */
template< class Scalar, class GridGeometry >
class CvfeFluxVariablesCache
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GridGeometry::FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
    using JacobianInverseTransposed = typename Element::Geometry::JacobianInverseTransposed;

public:

    //! update the cache for an scvf
    template< class Problem, class ElementVolumeVariables >
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace& scvf)
    {
        const auto geometry = element.geometry();
        const auto& localBasis = fvGeometry.feLocalBasis();

        // evaluate shape functions and gradients at the integration point
        ipGlobal_ = scvf.ipGlobal();
        const auto ipLocal = geometry.local(ipGlobal_);
        jacInvT_ = geometry.jacobianInverseTransposed(ipLocal);
        localBasis.evaluateJacobian(ipLocal, shapeJacobian_);
        localBasis.evaluateFunction(ipLocal, shapeValues_); // shape values for rho

        std::vector<ShapeValue> shapeValuesC;
        localBasis.evaluateFunction(geometry.local(element.geometry().center()), shapeValuesC); // shape values for rho

        // We add the contribution of the bubble function
        ShapeValue bubbleVal(1.0);
        for(auto v : shapeValues_)
            bubbleVal *= v;

        ShapeValue bubbleValC(1.0);
        for(auto v : shapeValuesC)
            bubbleValC *= v;

        for(int i=0; i<shapeValues_.size(); i++)
            shapeValues_[i] -= shapeValuesC[i]*bubbleVal/bubbleValC;

        ShapeJacobian bubbleJac(0.0);
        for(int i=0; i<shapeJacobian_.size(); i++)
        {
            ShapeJacobian jacVal = shapeJacobian_[i];
            for(int j=0; j<shapeValues_.size(); j++)
                if( i!= j)
                    jacVal *= shapeValues_[j];
            bubbleJac += jacVal;
        }
        bubbleJac /= bubbleValC;

        for(int i=0; i<shapeJacobian_.size(); i++)
        {
            ShapeJacobian val = bubbleJac;
            val *= shapeValuesC[i];
            shapeJacobian_[i] -= val;
        }

        shapeValues_.push_back(bubbleVal/bubbleValC);
        shapeJacobian_.push_back(bubbleJac);

        // compute the gradN at for every scv/dof
        gradN_.resize(fvGeometry.numScv());
        for (const auto& scv: scvs(fvGeometry))
            jacInvT_.mv(shapeJacobian_[scv.localDofIndex()][0], gradN_[scv.indexInElement()]);
    }

    //! returns the global position for which this cache has been updated
    const GlobalPosition& ipGlobal() const { return ipGlobal_; }
    //! returns the shape function gradients in local coordinates at the integration point
    const std::vector<ShapeJacobian>& shapeJacobian() const { return shapeJacobian_; }
    //! returns the shape function values at the integration point
    const std::vector<ShapeValue>& shapeValues() const { return shapeValues_; }
    //! returns inverse transposed jacobian at the integration point
    const JacobianInverseTransposed& jacInvT() const { return jacInvT_; }
    //! returns the shape function gradients in global coordinates at the local dof location
    const GlobalPosition& gradN(unsigned int localDofIndex) const { return gradN_[localDofIndex]; }

private:
    GlobalPosition ipGlobal_;
    std::vector<GlobalPosition> gradN_;
    std::vector<ShapeJacobian> shapeJacobian_;
    std::vector<ShapeValue> shapeValues_;
    JacobianInverseTransposed jacInvT_;
};

} // end namespace Dumux

#endif
