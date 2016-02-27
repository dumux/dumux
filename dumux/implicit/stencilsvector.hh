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
#ifndef DUMUX_STENCILSVECTOR_HH
#define DUMUX_STENCILSVECTOR_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{
/*!
 * \ingroup ImplicitModel
 * \brief Base class for a sub control volume, i.e a part of the control
 *        volume we are making the balance for.
 */
template<class TypeTag>
class StencilsVector : public std::vector<typename GET_PROP_TYPE(TypeTag, Stencils)>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Stencils = typename GET_PROP_TYPE(TypeTag, Stencils);
public:
    void update(const Problem& problem)
    {
        this->resize(problem.gridView().size(0));
        for (const auto& element : elements(problem.gridView()))
        {
            auto eIdx = problem.elementMapper().index(element);
            (*this)[eIdx].update(problem, element);
        }
    }
};

} // end namespace

#endif
