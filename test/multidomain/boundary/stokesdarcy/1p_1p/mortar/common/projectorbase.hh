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
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
#ifndef DUMUX_MORTAR_PROJECTOR_BASE_HH
#define DUMUX_MORTAR_PROJECTOR_BASE_HH

#include <string>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/promotiontraits.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/geometry/boundingboxtree.hh>
#include <dumux/common/geometry/geometricentityset.hh>

#include <dumux/discretization/method.hh>
#include <dumux/multidomain/embedded/mixeddimensionglue.hh>

namespace Dumux {

// namespace with some implementation details
namespace MortarProjectorDetail {

// Container to store pairs of dofIndices/shape values
template<class GridGeometry, class Scalar>
struct SubDomainShapeFunctionValues
{
    using IndexType = typename IndexTraits<typename GridGeometry::GridView>::GridIndex;
    using DofToValuePair = std::pair< IndexType, Scalar >;
    using Values = std::vector< DofToValuePair >;

    Values values;
};

// evaluates the shape values within an element at a given position (overload for box)
template<class Scalar, class GridGeometry,
         std::enable_if_t<GridGeometry::discMethod == DiscretizationMethod::box, int> = 0 >
SubDomainShapeFunctionValues<GridGeometry, Scalar>
evalShapeValues(const GridGeometry& gridGeometry,
                const typename GridGeometry::GridView::template Codim<0>::Entity& element,
                const typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate& globalPos)
{
    static constexpr int dim = GridGeometry::GridView::dimension;

    std::vector< Dune::FieldVector<double,1> > shapeValues;
    const auto& finiteElement = gridGeometry.feCache().get(element.geometry().type());
    const auto ipLocal = element.geometry().local( globalPos );
    finiteElement.localBasis().evaluateFunction(ipLocal, shapeValues);

    using ResultType = SubDomainShapeFunctionValues<GridGeometry, Scalar>;
    ResultType result;
    result.values.resize(finiteElement.localBasis().size());

    for (size_t i = 0; i < finiteElement.localBasis().size(); i++)
    {
        const auto dofIdx = gridGeometry.vertexMapper().subIndex(element, i, dim);
        result.values[i] = typename ResultType::DofToValuePair(dofIdx, shapeValues[i]);
    }

    return result;
}
} // end namespace MortarProjectorDetail

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template< class MortarFEBasis, class MortarSolutionVector,
          class SubDomainGridGeometry, class SubDomainSolutionVector >
class MortarProjectorBase
{
    using MortarGridView = typename MortarFEBasis::GridView;
    using MortarElement = typename MortarGridView::template Codim<0>::Entity;
    using MortarElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<MortarGridView>;
    using MortarGridIndex = typename MortarElementMapper::Index;
    using MortarScalar = typename MortarSolutionVector::field_type;

    using SubDomainGridView = typename SubDomainGridGeometry::GridView;
    using SubDomainElement = typename SubDomainGridView::template Codim<0>::Entity;
    using SubDomainElementMapper = typename SubDomainGridGeometry::ElementMapper;
    using SubDomainGridIndex = typename IndexTraits<SubDomainGridView>::GridIndex;
    using SubDomainScalar = typename SubDomainSolutionVector::field_type;

    static constexpr int mortarDim = MortarGridView::dimension;
    static constexpr int subDomainDim = SubDomainGridView::dimension;
    static_assert(mortarDim == subDomainDim - 1, "It must be mortarDim = sumDomainDum - 1!");

    using ScalarType = typename Dune::PromotionTraits<MortarScalar, SubDomainScalar>::PromotedType;

    // Structure containing data required for integration of a variable living
    // on one of the grids as seen from an element of the other grid
    template<class OtherGridIndex, class OtherScalar>
    struct IntegrationData
    {
        ScalarType overlapArea = 0;
        std::unordered_map<OtherGridIndex, OtherScalar> otherDofToWeightMap;
        std::vector<OtherGridIndex> coupledElementIndices;
    };

public:
    //! export solution vector of mortar domain
    using MortarSolution = MortarSolutionVector;
    //! export type used for scalar values
    using Scalar = ScalarType;

    //! The constructor
    MortarProjectorBase(std::shared_ptr<const MortarFEBasis> mortarFEBasis,
                        std::shared_ptr<const SubDomainGridGeometry> subDomainGridGeometry,
                        const std::string& paramGroup = "")
    : mortarFEBasis_(mortarFEBasis)
    , subDomainGridGeometry_(subDomainGridGeometry)
    , mortarElementMapper_(mortarFEBasis->gridView(), Dune::mcmgElementLayout())
    {
        // bounding box tree for mortar grid
        using MortarElementSet = GridViewGeometricEntitySet<MortarGridView, 0, MortarElementMapper>;
        using MortarBoundingBoxTree = Dumux::BoundingBoxTree<MortarElementSet>;
        MortarBoundingBoxTree mortarBBoxTree( std::make_shared<MortarElementSet>(mortarFEBasis->gridView(), mortarElementMapper_) );

        // intersect the bounding box trees
        using GlueType = MixedDimensionGlue< SubDomainGridView, MortarGridView, SubDomainElementMapper, MortarElementMapper >;
        GlueType glue(subDomainGridGeometry->boundingBoxTree(), mortarBBoxTree);

        std::size_t isCounter = 0;
        for (const auto& is : intersections(glue))
        {
            // we expect to have found only one neighbor per subdomain
            if (is.neighbor(0) != 1)
                DUNE_THROW(Dune::InvalidStateException, "Found " << is.neighbor(0) << " instead of one neighbor on mortar side");
            if (is.neighbor(1) != 1)
                DUNE_THROW(Dune::InvalidStateException, "Found " << is.neighbor(1) << " instead of one neighbor on sub-domain side");

            isCounter++;
            // "inside" element is mortar element
            const auto& mortarElement = is.inside(0);
            const auto& mortarElementGeometry = mortarElement.geometry();
            const auto mortarElementIdx = mortarElementMapper_.index(mortarElement);

            // "outside" is sub-domain element
            const auto& subDomainElement = is.outside(0);
            const auto subDomainElementIdx = subDomainGridGeometry->elementMapper().index(subDomainElement);

            // prepare entries for projections
            auto& mortarToSubDomainEntry = mortarToSubDomainProjection_[subDomainElementIdx];
            auto& subDomainToMortarEntry = subDomainToMortarProjection_[mortarElementIdx];

            const auto isGeometry = is.geometry();
            mortarToSubDomainEntry.overlapArea += isGeometry.volume();
            subDomainToMortarEntry.overlapArea += isGeometry.volume();

            // bind local view of mortar element
            auto localView = mortarFEBasis->localView();
            localView.bind(mortarElement);

            // perform integration over intersection geometry
            const auto order = getParamFromGroup<Scalar>(paramGroup, "ProjectionIntegrationOrder");
            const auto& quad = Dune::QuadratureRules<Scalar, mortarDim>::rule(isGeometry.type(), order);
            for (auto&& qp : quad)
            {
                const auto weight = qp.weight();
                const auto ie = isGeometry.integrationElement(qp.position());

                // Evaluate mortar shape function values at this point
                const auto& globalPos = isGeometry.global(qp.position());
                const auto ipLocalMortar = mortarElementGeometry.local( globalPos );
                const auto& mortarFiniteElement = localView.tree().finiteElement();

                std::vector< Dune::FieldVector<double,1> > mortarShapeFunctionValues;
                mortarFiniteElement.localBasis().evaluateFunction(ipLocalMortar, mortarShapeFunctionValues);

                // Evaluate sub-domain shape function values at this point
                const auto subDomainShapeValues = MortarProjectorDetail::evalShapeValues<Scalar>(*subDomainGridGeometry, subDomainElement, globalPos);

                // lambda to check if dofIdx is contained in map
                auto containsIdx = [] (const auto& map, const auto idx)
                {
                    auto it = std::find_if(map.begin(),
                                           map.end(),
                                           [idx] (const auto& pair) { return pair.first == idx; });
                    return it != map.end();
                };

                // store integration data
                for (size_t i = 0; i < mortarFiniteElement.localBasis().size(); i++)
                {
                    const auto mortarDofIdx = localView.index(i);

                    // if entry for dof doesn't exist yet, initialize to zero
                    if ( !containsIdx(mortarToSubDomainEntry.otherDofToWeightMap, mortarDofIdx) )
                        mortarToSubDomainEntry.otherDofToWeightMap[mortarDofIdx] = 0.0;
                    mortarToSubDomainEntry.otherDofToWeightMap[mortarDofIdx] += weight*ie*mortarShapeFunctionValues[i];

                    if ( !std::count(mortarToSubDomainEntry.coupledElementIndices.begin(),
                                     mortarToSubDomainEntry.coupledElementIndices.end(),
                                     mortarElementIdx) )
                        mortarToSubDomainEntry.coupledElementIndices.push_back(mortarElementIdx);
                }

                for (size_t i = 0; i < subDomainShapeValues.values.size(); i++)
                {
                    const auto subDomainDofIdx = subDomainShapeValues.values[i].first;

                    // if entry for dof doesn't exist yet, initialize to zero
                    if ( !containsIdx(subDomainToMortarEntry.otherDofToWeightMap, subDomainDofIdx) )
                        subDomainToMortarEntry.otherDofToWeightMap[subDomainDofIdx] = 0.0;
                    subDomainToMortarEntry.otherDofToWeightMap[subDomainDofIdx] += weight*ie*subDomainShapeValues.values[i].second;

                    if ( !std::count(subDomainToMortarEntry.coupledElementIndices.begin(),
                                     subDomainToMortarEntry.coupledElementIndices.end(),
                                     subDomainElementIdx) )
                        subDomainToMortarEntry.coupledElementIndices.push_back(subDomainElementIdx);
                }
            }
        }

        std::cout << isCounter << " intersections were found between mortar and sub-domain grid" << std::endl;
    }

    //! set mortar sub-domain solution vector
    void setMortarSolutionPointer(std::shared_ptr<const MortarSolutionVector> x)
    {
        mortarSolution_ = x;
    }

    //! set sub-domain solution vector
    void setSubDomainSolutionPointer(std::shared_ptr<const SubDomainSolutionVector> x)
    {
        subDomainSolution_ = x;
    }

    //! returns true if the element has an overlap with the mortar grid
    bool hasOverlapWithMortar(const SubDomainElement& subDomainElement) const
    {
        const auto subDomainElementIdx = subDomainGridGeometry_->elementMapper().index(subDomainElement);
        return mortarToSubDomainProjection_.find(subDomainElementIdx) != mortarToSubDomainProjection_.end();
    }

    //! integrates the mortar variable over the overlap with a sub-domain element
    MortarScalar integrateMortarVariable(const SubDomainElement& subDomainElement) const
    {
        assert(hasOverlapWithMortar(subDomainElement));
        const auto subDomainElementIdx = subDomainGridGeometry_->elementMapper().index(subDomainElement);
        const auto& projectionEntry = mortarToSubDomainProjection_.at(subDomainElementIdx);

        MortarScalar value = 0.0;
        for (const auto& entry : projectionEntry.otherDofToWeightMap)
            value += (*mortarSolution_)[entry.first]*entry.second;

        return value;
    }

    //! returns the overlap area between sub-domain and mortar domain
    Scalar overlapArea(const SubDomainElement& subDomainElement) const
    {
        assert(hasOverlapWithMortar(subDomainElement));
        const auto subDomainElementIdx = subDomainGridGeometry_->elementMapper().index(subDomainElement);
        const auto& projectionEntry = mortarToSubDomainProjection_.at(subDomainElementIdx);
        return projectionEntry.overlapArea;
    }

protected:
    std::shared_ptr<const MortarFEBasis> mortarFEBasis_;
    std::shared_ptr<const SubDomainGridGeometry> subDomainGridGeometry_;

    std::shared_ptr<const MortarSolutionVector> mortarSolution_;
    std::shared_ptr<const SubDomainSolutionVector> subDomainSolution_;

    MortarElementMapper mortarElementMapper_;

    std::unordered_map< SubDomainGridIndex, IntegrationData<MortarGridIndex, MortarScalar> > mortarToSubDomainProjection_;
    std::unordered_map< MortarGridIndex, IntegrationData<SubDomainGridIndex, SubDomainScalar> > subDomainToMortarProjection_;
};

} // end namespace Dumux

#endif
