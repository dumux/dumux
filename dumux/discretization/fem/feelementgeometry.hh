// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FEMDiscretization
 * \brief Grid geometry local view, which is a wrapper around a
 *        finite element basis local view.
 */
#ifndef DUMUX_DISCRETIZATION_FE_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_FE_ELEMENT_GEOMETRY_HH

#include <utility>

namespace Dumux {

/*!
 * \ingroup FEMDiscretization
 * \brief Grid geometry local view, which is a wrapper around a
 *        finite element basis local view.
 * \tparam The grid geometry type
 * \tparam The FEBasis local view
 */
template<class GridGeometry>
class FEElementGeometry
{
    using GridView = typename GridGeometry::GridView;

    using FEBasis = typename GridGeometry::FEBasis;
    using FEBasisLocalView = typename FEBasis::LocalView;

public:
    //! export type of the element
    using Element = typename GridView::template Codim<0>::Entity;

    //! constructor taking grid geometry
    FEElementGeometry(const GridGeometry& gg)
    : gridGeometry_(gg)
    , feBasisLocalView_(gg.feBasis().localView())
    {}

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    FEElementGeometry bind(const Element& element) &&
    {
        this->bind(element);
        return std::move(*this);
    }

    //! Prepare element-local data
    void bind(const Element& element) &
    { feBasisLocalView_.bind(element); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    FEElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! Prepare element-local data
    void bindElement(const Element& element) &
    { bind(element); }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return feBasisLocalView_.isBound(); }

    //! The bound element
    const Element& element() const
    { return feBasisLocalView_.element(); }

    //! Return the finite element basis local view
    const FEBasisLocalView& feBasisLocalView() const
    { return feBasisLocalView_; }

    //! Return reference to the grid geometry
    const GridGeometry& gridGeometry() const
    { return gridGeometry_; }

private:
    const GridGeometry& gridGeometry_;
    FEBasisLocalView feBasisLocalView_;
};

} // end namespace Dumux

#endif
