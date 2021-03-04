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
 * \ingroup Discretization
 * \brief Class that contains the element-local (or element stencil-local) data
 *        required to evaluate the terms of discrete equations.
 */
#ifndef DUMUX_LOCAL_CONTEXT_HH
#define DUMUX_LOCAL_CONTEXT_HH

#include <type_traits>

namespace Dumux {
namespace Experimental {

class EmptyCouplingContext {};

/*!
 * \ingroup Discretization
 * \brief TODO: Doc me
 */
template<class EV, class CC = EmptyCouplingContext>
class LocalContext
{

    using CouplingContext = CC;
    using GridVariables = typename EV::GridVariables;
    using GridGeometry  = typename GridVariables::GridGeometry;
    using GridView = typename GridGeometry::GridView;

    static constexpr bool isEmptyCC = std::is_same_v<CC, EmptyCouplingContext>;

public:

    using Element = typename GridView::template Codim<0>::Entity;
    using ElementGridGeometry = typename GridGeometry::LocalView;
    using ElementVariables = EV;

    //! Constructor for contexts without coupling data
    template<bool e = isEmptyCC, std::enable_if_t<e, int> = 0>
    LocalContext(const Element& element,
                 const ElementGridGeometry& eg,
                 const ElementVariables& ev)
    : element_(&element)
    , elemGeom_(&eg)
    , elemVars_(&ev)
    {}

    //! Constructor for contexts with coupling data
    template<bool e = isEmptyCC, std::enable_if_t<!e, int> = 0>
    LocalContext(const Element& element,
                 const ElementGridGeometry& eg,
                 const ElementVariables& ev,
                 const CouplingContext& cc)
    : element_(&element)
    , elemGeom_(&eg)
    , elemVars_(&ev)
    , couplingContext_(&cc)
    {}

    //! TODO: Doc me
    const Element& element() const { return *element_; }
    const ElementGridGeometry elementGridGeometry() const { return *elemGeom_; }
    const ElementVariables& elementVariables() const { return *elemVars_; }

    //! TODO: Doc me
    template<bool e = isEmptyCC, std::enable_if_t<!e, int> = 0>
    const CouplingContext& couplingContext() const { return *couplingContext_; }

private:
    const Element* element_;
    const ElementGridGeometry* elemGeom_;
    const ElementVariables* elemVars_;
    const CouplingContext* couplingContext_;
};

/*!
 * \ingroup Discretization
 * \brief TODO: Doc me
 */
template<class EV>
LocalContext<EV>
makeLocalContext(const typename EV::GridVariables::GridGeometry::GridView::template Codim<0>::Entity& element,
                 const typename EV::GridVariables::GridGeometry::LocalView& gglocalView,
                 const EV& elemVariables)
{ return {element, gglocalView, elemVariables}; }

/*!
 * \ingroup Discretization
 * \brief TODO: Doc me
 */
template<class EV, class CC>
LocalContext<EV, CC>
makeLocalContext(const typename EV::GridVariables::GridGeometry::GridView::template Codim<0>::Entity& element,
                 const typename EV::GridVariables::GridGeometry::LocalView& gglocalView,
                 const EV& elemVariables,
                 const CC& couplingContext)
{
    static_assert(!std::is_same_v<CC, EmptyCouplingContext>, "Invalid coupling context!");
    return {element, gglocalView, elemVariables, couplingContext};
}

} // end namespace Experimental
} // end namespace Dumux

#endif
