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
 * \brief Implements the notion of stencils
 */
#ifndef DUMUX_CC_STENCILS_HH
#define DUMUX_CC_STENCILS_HH

#include <set>
#include <dumux/implicit/cellcentered/properties.hh>

namespace Dumux
{
/*!
 * \ingroup ImplicitModel
 * \brief Base class for a sub control volume, i.e a part of the control
 *        volume we are making the balance for.
 */
template<class TypeTag>
class CCStencils
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    // TODO a separate stencil class all stencils can derive from?
    using Stencil = std::set<IndexType>;
public:
    void update(const Problem& problem, const Element& element)
    {
        for (auto&& scvf : problem.model().fvGeometries(element).scvfs())
        {
            auto&& fluxStencil = problem.model().fluxVars(scvf).stencil();
            elementStencil_.insert(fluxStencil.begin(), fluxStencil.end());
        }
        neighborStencil_ = elementStencil_;
        neighborStencil_.erase(problem.elementMapper().index(element));
    }

    //! The full element stencil (all element this element is interacting with)
    const Stencil& elementStencil() const
    {
        return elementStencil_;
    }

    //! //! The full element stencil without this element
    const Stencil& neighborStencil() const
    {
        return neighborStencil_;
    }

private:
    Stencil elementStencil_;
    Stencil neighborStencil_;
};

} // end namespace

#endif
