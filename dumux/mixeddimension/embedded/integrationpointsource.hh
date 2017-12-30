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
 * \ingroup EmbeddedCoupling
 * \brief An integration point source class,
 *        i.e. sources located at a single point in space
 *        associated with a quadrature point
 */

#ifndef DUMUX_INTEGRATION_POINTSOURCE_HH
#define DUMUX_INTEGRATION_POINTSOURCE_HH

#include <dumux/common/pointsource.hh>

namespace Dumux
{

/*!
 * \ingroup EmbeddedCoupling
 * \brief An integration point source class with an identifier to attach data
 *        and a quadrature weight and integration element
 */
template<class TypeTag, typename IdType = std::size_t>
class IntegrationPointSource : public Dumux::IdPointSource<TypeTag, IdType>
{
    using ParentType = Dumux::IdPointSource<TypeTag, IdType>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    static const int dimworld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    //! Constructor for integration point sources
    IntegrationPointSource(GlobalPosition pos, PrimaryVariables values, IdType id,
                           Scalar qpweight, Scalar integrationElement,
                           const std::vector<std::size_t>& elementIndices)
      : ParentType(pos, values, id),
        qpweight_(qpweight), integrationElement_(integrationElement),
        elementIndices_(elementIndices) {}

    //! Constructor for integration point sources, when there is no
    // value known at the time of initialization
    IntegrationPointSource(GlobalPosition pos, IdType id,
                           Scalar qpweight, Scalar integrationElement,
                           const std::vector<std::size_t>& elementIndices)
      : ParentType(pos, id),
        qpweight_(qpweight), integrationElement_(integrationElement),
        elementIndices_(elementIndices) {}


    Scalar quadratureWeight() const
    {
        return qpweight_;
    }

    Scalar integrationElement() const
    {
        return integrationElement_;
    }

    const std::vector<std::size_t>& elementIndices() const
    {
        return elementIndices_;
    }

    //! Convenience = operator overload modifying only the values
    IntegrationPointSource& operator= (const PrimaryVariables& values)
    {
        ParentType::operator=(values);
        return *this;
    }

    //! Convenience = operator overload modifying only the values
    IntegrationPointSource& operator= (Scalar s)
    {
        ParentType::operator=(s);
        return *this;
    }

private:
    Scalar qpweight_;
    Scalar integrationElement_;
    std::vector<std::size_t> elementIndices_;
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief A helper class calculating a DOF-index to point source map
 */
class IntegrationPointSourceHelper
{

public:
    //! calculate a DOF index to point source map from given vector of point sources
    template<class FVGridGeometry, class PointSource>
    static void computePointSourceMap(const FVGridGeometry& fvGridGeometry,
                                      std::vector<PointSource>& sources,
                                      std::map<std::pair<unsigned int, unsigned int>, std::vector<PointSource> >& pointSourceMap)
    {
        for (auto&& source : sources)
        {
            // compute in which elements the point source falls
            const auto& entities = source.elementIndices();
            // split the source values equally among all concerned entities
            source.setEmbeddings(source.embeddings()*entities.size());
            // loop over all intersected elements
            for (unsigned int eIdx : entities)
            {
                // if(isBox)
                // {
                //     // check in which subcontrolvolume(s) we are
                //     const auto element = boundingBoxTree.entitySet().entity(eIdx);
                //     auto fvGeometry = localView(fvGridGeometry);
                //     fvGeometry.bindElement(element);
                //     const auto globalPos = source.position();
                //     // loop over all sub control volumes and check if the point source is inside
                //     std::vector<unsigned int> scvIndices;
                //     for (auto&& scv : scvs(fvGeometry))
                //     {
                //         if (intersectsPointGeometry(globalPos, scv.geometry()))
                //             scvIndices.push_back(scv.indexInElement());
                //     }
                //     // for all scvs that where tested positiv add the point sources
                //     // to the element/scv to point source map
                //     for (auto scvIdx : scvIndices)
                //     {
                //         const auto key = std::make_pair(eIdx, scvIdx);
                //         if (pointSourceMap.count(key))
                //             pointSourceMap.at(key).push_back(source);
                //         else
                //             pointSourceMap.insert({key, {source}});
                //         // split equally on the number of matched scvs
                //         auto& s = pointSourceMap.at(key).back();
                //         s.setEmbeddings(scvIndices.size()*s.embeddings());
                //     }
                // }
                // else
                // {
                    // add the pointsource to the DOF map
                    const auto key = std::make_pair(eIdx, /*scvIdx=*/ 0);
                    if (pointSourceMap.count(key))
                        pointSourceMap.at(key).push_back(source);
                    else
                        pointSourceMap.insert({key, {source}});
                // }
            }
        }
    }
};

} // end namespace Dumux

#endif
