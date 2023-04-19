// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief An integration point source class,
 *        i.e. sources located at a single point in space
 *        associated with a quadrature point
 */

#ifndef DUMUX_INTEGRATION_POINTSOURCE_HH
#define DUMUX_INTEGRATION_POINTSOURCE_HH

#include <type_traits>
#include <dune/common/reservedvector.hh>
#include <dumux/common/pointsource.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief An integration point source class with an identifier to attach data
 *        and a quadrature weight and integration element
 */
template<class GlobalPosition, class SourceValues, typename IdType = std::size_t>
class IntegrationPointSource : public IdPointSource<GlobalPosition, SourceValues, IdType>
{
    using ParentType = IdPointSource<GlobalPosition, SourceValues, IdType>;
    using Scalar = std::decay_t<decltype(std::declval<SourceValues>()[0])>;

public:
    //! Constructor for integration point sources
    IntegrationPointSource(GlobalPosition pos, SourceValues values, IdType id,
                           Scalar qpweight, Scalar integrationElement,
                           std::size_t elementIndex)
    : ParentType(pos, values, id)
    , qpweight_(qpweight)
    , integrationElement_(integrationElement)
    , elementIndex_(elementIndex)
    {}

    //! Constructor for integration point sources, when there is no
    //! value known at the time of initialization
    IntegrationPointSource(GlobalPosition pos, IdType id,
                           Scalar qpweight, Scalar integrationElement,
                           std::size_t elementIndex)
    : ParentType(pos, id)
    , qpweight_(qpweight)
    , integrationElement_(integrationElement)
    , elementIndex_(elementIndex)
    {}

    Scalar quadratureWeight() const
    { return qpweight_; }

    Scalar integrationElement() const
    { return integrationElement_; }

    void setQuadratureWeight(const Scalar qpWeight)
    { return qpweight_ = qpWeight; }

    void setIntegrationElement(const Scalar ie)
    { integrationElement_ = ie; }

    //! The index of the element this integration point source is associated with
    std::size_t elementIndex() const
    { return elementIndex_; }

    //! Convenience = operator overload modifying only the values
    IntegrationPointSource& operator= (const SourceValues& values)
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
    std::size_t elementIndex_;
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief A helper class calculating a DOF-index to point source map
 */
class IntegrationPointSourceHelper
{

public:
    //! calculate a DOF index to point source map from given vector of point sources
    template<class GridGeometry, class PointSource, class PointSourceMap>
    static void computePointSourceMap(const GridGeometry& gridGeometry,
                                      const std::vector<PointSource>& sources,
                                      PointSourceMap& pointSourceMap,
                                      const std::string& paramGroup = "")
    {
        for (const auto& source : sources)
        {
            // get the index of the element in which the point source falls
            const auto eIdx = source.elementIndex();
            if constexpr (GridGeometry::discMethod == DiscretizationMethods::box
                || GridGeometry::discMethod == DiscretizationMethods::fcdiamond)
            {
                auto fvGeometry = localView(gridGeometry);
                // check in which subcontrolvolume(s) we are
                const auto element = gridGeometry.boundingBoxTree().entitySet().entity(eIdx);
                fvGeometry.bindElement(element);
                const auto globalPos = source.position();

                static const bool boxPointSourceLumping = getParamFromGroup<bool>(paramGroup, "PointSource.EnableBoxLumping", true);
                if (boxPointSourceLumping)
                {
                    // loop over all sub control volumes and check if the point source is inside
                    constexpr int dim = GridGeometry::GridView::dimension;
                    Dune::ReservedVector<std::size_t, 1<<dim> scvIndices;
                    for (const auto& scv : scvs(fvGeometry))
                        if (intersectsPointGeometry(globalPos, fvGeometry.geometry(scv)))
                            scvIndices.push_back(scv.indexInElement());

                    // for all scvs that where tested positive add the point sources
                    // to the element/scv to point source map
                    for (auto scvIdx : scvIndices)
                    {
                        const auto key = std::make_pair(eIdx, scvIdx);
                        if (pointSourceMap.count(key))
                            pointSourceMap.at(key).push_back(source);
                        else
                            pointSourceMap.insert({key, {source}});
                        // split equally on the number of matched scvs
                        auto& s = pointSourceMap.at(key).back();
                        s.setEmbeddings(scvIndices.size()*s.embeddings());
                    }
                }

                // distribute the sources using the basis function weights
                else
                {
                    const auto& localBasis = fvGeometry.feLocalBasis();
                    const auto ipLocal = element.geometry().local(globalPos);
                    using Scalar = std::decay_t<decltype(source.values()[0])>;
                    std::vector<typename Dune::FieldVector<Scalar, 1>> shapeValues;
                    localBasis.evaluateFunction(ipLocal, shapeValues);
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        const auto key = std::make_pair(eIdx, scv.indexInElement());
                        if (pointSourceMap.count(key))
                            pointSourceMap.at(key).push_back(source);
                        else
                            pointSourceMap.insert({key, {source}});

                        // adjust the integration element
                        auto& s = pointSourceMap.at(key).back();
                        s.setIntegrationElement(shapeValues[scv.indexInElement()]*s.integrationElement());
                    }
                }
            }
            else
            {
                // add the pointsource to the DOF map
                const auto key = std::make_pair(eIdx, /*scvIdx=*/ 0);
                if (pointSourceMap.count(key))
                    pointSourceMap.at(key).push_back(source);
                else
                    pointSourceMap.insert({key, {source}});
            }
        }
    }
};

} // end namespace Dumux

#endif
