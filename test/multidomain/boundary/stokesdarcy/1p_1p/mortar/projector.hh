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
#ifndef DUMUX_MORTAR_DARCY_PROJECTOR_HH
#define DUMUX_MORTAR_DARCY_PROJECTOR_HH

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

#include <dumux/multidomain/embedded/mixeddimensionglue.hh>

namespace Dumux {

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template< class MortarFEBasis, class MortarSolutionVector,
          class DarcyGridGeometry, class DarcySolutionVector >
class MortarDarcyProjector
{
    using MortarGridView = typename MortarFEBasis::GridView;
    using MortarElement = typename MortarGridView::template Codim<0>::Entity;
    using MortarElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<MortarGridView>;
    using MortarGridIndex = typename MortarElementMapper::Index;
    using MortarScalar = typename MortarSolutionVector::field_type;

    using DarcyGridView = typename DarcyGridGeometry::GridView;
    using DarcyElement = typename DarcyGridView::template Codim<0>::Entity;
    using DarcyElementMapper = typename DarcyGridGeometry::ElementMapper;
    using DarcyGridIndex = typename IndexTraits<DarcyGridView>::GridIndex;
    using DarcyScalar = typename DarcySolutionVector::field_type;

    static constexpr int mortarDim = MortarGridView::dimension;
    static constexpr int darcyDim = DarcyGridView::dimension;
    static_assert(mortarDim == darcyDim - 1, "It must be mortarDim = darcyDim - 1!");

    using Scalar = typename Dune::PromotionTraits<MortarScalar, DarcyScalar>::PromotedType;

    // Structure containing data required for integration of a variable living
    // on one of the grids as seen from an element of the other grid
    template<class OtherGridIndex, class OtherScalar>
    struct IntegrationData
    {
        Scalar overlapArea = 0;
        std::unordered_map<OtherGridIndex, OtherScalar> otherDofToWeightMap;
    };

public:
    //! export solution vector of mortar domain
    using MortarSolution = MortarSolutionVector;

    //! The constructor
    MortarDarcyProjector(std::shared_ptr<const MortarFEBasis> mortarFEBasis,
                         std::shared_ptr<const DarcyGridGeometry> darcyGridGeometry,
                         const std::string& paramGroup = "")
    : mortarFEBasis_(mortarFEBasis)
    , darcyGridGeometry_(darcyGridGeometry)
    , mortarElementMapper_(mortarFEBasis->gridView(), Dune::mcmgElementLayout())
    {
        // bounding box tree for mortar grid
        using MortarElementSet = GridViewGeometricEntitySet<MortarGridView, 0, MortarElementMapper>;
        using MortarBoundingBoxTree = Dumux::BoundingBoxTree<MortarElementSet>;
        MortarBoundingBoxTree mortarBBoxTree( std::make_shared<MortarElementSet>(mortarFEBasis->gridView(), mortarElementMapper_) );

        // intersect the bounding box trees
        using GlueType = MixedDimensionGlue< DarcyGridView, MortarGridView, DarcyElementMapper, MortarElementMapper >;
        GlueType glue(darcyGridGeometry->boundingBoxTree(), mortarBBoxTree);

        std::size_t isCounter = 0;
        for (const auto& is : intersections(glue))
        {
            // we expect to have found only one neighbor per subdomain
            if (is.neighbor(0) != 1)
                DUNE_THROW(Dune::InvalidStateException, "Found " << is.neighbor(0) << " instead of one neighbor on mortar side");
            if (is.neighbor(1) != 1)
                DUNE_THROW(Dune::InvalidStateException, "Found " << is.neighbor(1) << " instead of one neighbor on Darcy side");

            isCounter++;
            // "inside" element is mortar element
            const auto& mortarElement = is.inside(0);
            const auto& mortarElementGeometry = mortarElement.geometry();
            const auto mortarElementIdx = mortarElementMapper_.index(mortarElement);

            // "outside" is Darcy element
            const auto& darcyElement = is.outside(0);
            const auto& darcyElementGeometry = darcyElement.geometry();
            const auto darcyElementIdx = darcyGridGeometry->elementMapper().index(darcyElement);

            // prepare entries for projections
            auto& mortarToDarcyEntry = mortarToDarcyProjection_[darcyElementIdx];
            auto& darcyToMortarEntry = darcyToMortarProjection_[mortarElementIdx];

            const auto isGeometry = is.geometry();
            mortarToDarcyEntry.overlapArea += isGeometry.volume();
            darcyToMortarEntry.overlapArea += isGeometry.volume();

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

                // Evaluate mortar/darcy shape function values at this point
                std::vector< Dune::FieldVector<double,1> > mortarShapeFunctionValues;
                std::vector< Dune::FieldVector<double,1> > darcyShapeFunctionValues;

                const auto& mortarFiniteElement = localView.tree().finiteElement();
                const auto& darcyFiniteElement = darcyGridGeometry->feCache().get(darcyElement.geometry().type());

                const auto& globalPos = isGeometry.global(qp.position());
                const auto ipLocalMortar = mortarElementGeometry.local( globalPos );
                const auto ipLocalDarcy = darcyElementGeometry.local( globalPos );
                mortarFiniteElement.localBasis().evaluateFunction(ipLocalMortar, mortarShapeFunctionValues);
                darcyFiniteElement.localBasis().evaluateFunction(ipLocalDarcy, darcyShapeFunctionValues);

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
                    if ( !containsIdx(mortarToDarcyEntry.otherDofToWeightMap, mortarDofIdx) )
                        mortarToDarcyEntry.otherDofToWeightMap[mortarDofIdx] = 0.0;

                    mortarToDarcyEntry.otherDofToWeightMap[mortarDofIdx] += weight*ie*mortarShapeFunctionValues[i];
                }

                for (size_t i = 0; i < darcyFiniteElement.localBasis().size(); i++)
                {
                    const auto darcyDofIdx = darcyGridGeometry->vertexMapper().subIndex(darcyElement, i, darcyDim);

                    // if entry for dof doesn't exist yet, initialize to zero
                    if ( !containsIdx(darcyToMortarEntry.otherDofToWeightMap, darcyDofIdx) )
                        darcyToMortarEntry.otherDofToWeightMap[darcyDofIdx] = 0.0;

                    darcyToMortarEntry.otherDofToWeightMap[darcyDofIdx] += weight*ie*darcyShapeFunctionValues[i];
                }
            }
        }

        std::cout << isCounter << " intersections were found between mortar and darcy grid" << std::endl;
    }

    //! set mortar sub-domain solution vector
    void setMortarSolutionPointer(std::shared_ptr<const MortarSolutionVector> x)
    {
        mortarSolution_ = x;
    }

    //! set darcy sub-domain solution vector
    void setDarcySolutionPointer(std::shared_ptr<const DarcySolutionVector> x)
    {
        darcySolution_ = x;
    }

    //! returns true if the element has an overlap with the mortar grid
    bool hasOverlapWithMortar(const DarcyElement& darcyElement) const
    {
        const auto darcyElementIdx = darcyGridGeometry_->elementMapper().index(darcyElement);
        return mortarToDarcyProjection_.find(darcyElementIdx) != mortarToDarcyProjection_.end();
    }

    //! integrates the mortar variable over the overlap with a darcy element
    MortarScalar integrateMortarVariable(const DarcyElement& darcyElement) const
    {
        assert(hasOverlapWithMortar(darcyElement));
        const auto darcyElementIdx = darcyGridGeometry_->elementMapper().index(darcyElement);
        const auto& projectionEntry = mortarToDarcyProjection_.at(darcyElementIdx);

        MortarScalar value = 0.0;
        for (const auto& entry : projectionEntry.otherDofToWeightMap)
            value += (*mortarSolution_)[entry.first]*entry.second;

        return value;
    }

    //! returns the overlap area between darcy and mortar domain
    Scalar overlapArea(const DarcyElement& darcyElement) const
    {
        assert(hasOverlapWithMortar(darcyElement));
        const auto darcyElementIdx = darcyGridGeometry_->elementMapper().index(darcyElement);
        const auto& projectionEntry = mortarToDarcyProjection_.at(darcyElementIdx);
        return projectionEntry.overlapArea;
    }

    //! projects the interface pressures in the Darcy domain to mortar space
    MortarSolutionVector projectInterfacePressures() const
    {
        MortarSolutionVector p;
        p.resize(mortarFEBasis_->gridView().size(0));

        for (const auto& element : elements(mortarFEBasis_->gridView()))
        {
            const auto& eIdx = mortarElementMapper_.index(element);
            const auto& projectionEntry = darcyToMortarProjection_.at(eIdx);

            // integrate pressure acting from darcy domain on this element
            DarcyScalar elemPressure = 0.0;
            for (const auto& entry : projectionEntry.otherDofToWeightMap)
                elemPressure += (*darcySolution_)[entry.first]*entry.second;
            p[eIdx] = elemPressure;
        }

        return p;
    }

private:
    std::shared_ptr<const MortarFEBasis> mortarFEBasis_;
    std::shared_ptr<const DarcyGridGeometry> darcyGridGeometry_;

    std::shared_ptr<const MortarSolutionVector> mortarSolution_;
    std::shared_ptr<const DarcySolutionVector> darcySolution_;

    MortarElementMapper mortarElementMapper_;

    std::unordered_map< DarcyGridIndex, IntegrationData<MortarGridIndex, MortarScalar> > mortarToDarcyProjection_;
    std::unordered_map< MortarGridIndex, IntegrationData<DarcyGridIndex, DarcyScalar> > darcyToMortarProjection_;
};

} // end namespace Dumux

#endif
