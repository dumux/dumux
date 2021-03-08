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

namespace Dumux {
namespace Experimental {

/*!
 * \ingroup Discretization
 * \brief Implementation of element-stencil-local contexts, which solely store
 *        the local geometry and primary/secondary variables. This implementation
 *        defines the minimum interface required for contexts in single-domain
 *        applications to work with the default assembly mechanism.
 * \tparam EV The element-local view on the grid variables
 */
template<class EV>
class DefaultLocalContext
{
    using GridVariables = typename EV::GridVariables;
    using GridGeometry  = typename GridVariables::GridGeometry;

public:
    using ElementGridGeometry = typename GridGeometry::LocalView;
    using ElementVariables = EV;

    //! Constructor
    DefaultLocalContext(const ElementGridGeometry& eg,
                        const ElementVariables& ev)
    : elemGeom_(&eg)
    , elemVars_(&ev)
    {}

    //! Return the element-local view on the grid geometry
    const ElementGridGeometry elementGridGeometry() const
    { return *elemGeom_; }

    //! Return the element-local view on the grid variables
    const ElementVariables& elementVariables() const
    { return *elemVars_; }

private:
    const ElementGridGeometry* elemGeom_;
    const ElementVariables* elemVars_;
};

/*!
 * \ingroup Discretization
 * \brief Factory function to create a context from local views.
 */
template<class EV>
DefaultLocalContext<EV>
makeLocalContext(const typename EV::GridVariables::GridGeometry::LocalView& gglocalView,
                 const EV& elemVariables)
{ return {gglocalView, elemVariables}; }

} // end namespace Experimental
} // end namespace Dumux

#endif
