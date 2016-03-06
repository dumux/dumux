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
 * \brief Base class for the volume variables vector
 */
#ifndef DUMUX_IMPLICIT_VOLVARSVECTOR_HH
#define DUMUX_IMPLICIT_VOLVARSVECTOR_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the volume variables vector
 */
template<class TypeTag>
class VolumeVariablesVector : public std::vector<typename GET_PROP_TYPE(TypeTag, VolumeVariables)>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);

public:
    void update(const Problem& problem, const SolutionVector& sol)
    {
        for (const auto& element : elements(problem.gridView()))
        {
            for (auto&& scv : problem.model().fvGeometries(element).scvs())
            {
                (*this)[scv.index()].update(sol[scv.dofIndex()],
                                            problem,
                                            element,
                                            scv);
            }

        }
    }
};

} // end namespace

#endif
