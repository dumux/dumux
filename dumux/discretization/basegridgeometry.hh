// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Base class for grid geometries
 */
#ifndef DUMUX_DISCRETIZATION_BASE_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_BASE_GRID_GEOMETRY_HH

#include <memory>
#include <dune/common/std/type_traits.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/basicgridgeometry.hh>

// make the local view function available whenever we use the grid geometry
#include <dumux/discretization/localview.hh>

namespace Dumux {

namespace Detail {
template<class T>
using SpecifiesBaseGridGeometry = typename T::BasicGridGeometry;

template<class T>
using SpecifiesGeometryHelper = typename T::GeometryHelper;
} // end namespace Detail

/*!
 * \ingroup Discretization
 * \brief Type of the basic grid geometry implementation used as backend
 */
template<class GV, class T>
using BasicGridGeometry_t = Dune::Std::detected_or_t<
    Dumux::BasicGridGeometry<GV, typename T::ElementMapper, typename T::VertexMapper>,
    Detail::SpecifiesBaseGridGeometry,
    T
>;

/*!
 * \ingroup Discretization
 * \brief Base class for all grid geometries
 * \tparam GV the grid view type
 * \tparam Traits traits class that specifies mappers and basic grid geometry
 */
template<class GV, class Traits>
class BaseGridGeometry
{
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using Element = typename GV::template Codim<0>::Entity;
    using BaseImplementation = BasicGridGeometry_t<GV, Traits>;
public:
    //! export the grid type
    using Grid = typename BaseImplementation::Grid;
    //! export the grid view type
    using GridView = typename BaseImplementation::GridView;
    //! export the global coordinate type
    using GlobalCoordinate = typename BaseImplementation::GlobalCoordinate;
    //! export the element mapper type
    using ElementMapper = typename BaseImplementation::ElementMapper;
    //! export the vertex mapper type
    using VertexMapper = typename BaseImplementation::VertexMapper;

    /*!
     * \ingroup Discretization
     * \brief Constructor from a BaseImplementation
     */
    BaseGridGeometry(std::shared_ptr<BaseImplementation> impl)
    : impl_(std::move(impl))
    {}

    /*!
     * \ingroup Discretization
     * \brief Constructor from a grid view
     * \param gridView the grid view on which to construct the grid geometry
     */
    BaseGridGeometry(const GridView& gridView)
    : BaseGridGeometry(std::make_shared<BaseImplementation>(gridView))
    {}

    /*!
     * \brief Update all fvElementGeometries (call this after grid adaption)
     */
    void update(const GridView& gridView)
    { impl_->update(gridView); }

    /*!
     * \brief Update all fvElementGeometries (call this after grid adaption)
     */
    void update(GridView&& gridView)
    { impl_->update(std::move(gridView)); }

    /*!
     * \brief Return the gridView this grid geometry object lives on
     */
    const GridView& gridView() const
    { return impl_->gridView(); }

    /*!
     * \brief Returns the mapper for vertices to indices for constant grids.
     */
    const VertexMapper &vertexMapper() const
    { return impl_->vertexMapper(); }

    /*!
     * \brief Returns the mapper for elements to indices for constant grids.
     */
    const ElementMapper &elementMapper() const
    { return impl_->elementMapper(); }

    /*!
     * \brief Returns the mapper for vertices to indices for possibly adaptive grids.
     */
    VertexMapper &vertexMapper()
    { return impl_->vertexMapper(); }

    /*!
     * \brief Returns the mapper for elements to indices for possibly adaptive grids.
     */
    ElementMapper &elementMapper()
    { return impl_->elementMapper(); }

    /*!
     * \brief Returns the bounding box tree of the grid
     */
    decltype(auto) boundingBoxTree() const
    { return impl_->boundingBoxTree(); }

    /*!
     * \brief Returns the element index to element map
     */
    decltype(auto) elementMap() const
    { return impl_->elementMap(); }

    /*!
     * \brief Get an element from a global element index
     */
    Element element(GridIndexType eIdx) const
    { return impl_->element(eIdx); }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the smallest values.
     */
    const GlobalCoordinate &bBoxMin() const
    { return impl_->bBoxMin(); }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the largest values.
     */
    const GlobalCoordinate &bBoxMax() const
    { return impl_->bBoxMax(); }

    /*!
     * \brief Returns if the grid geometry is periodic (at all)
     */
    bool isPeriodic() const
    { return periodic_; }

protected:
    /*!
     * \brief Set the periodicity of the grid geometry
     */
    void setPeriodic(bool value = true)
    { periodic_ = value; }

private:
    std::shared_ptr<BaseImplementation> impl_;

    //! if the grid geometry has periodic boundaries
    bool periodic_ = false;
};

} // end namespace Dumux

#endif
