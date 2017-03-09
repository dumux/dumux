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
 * \brief Base class for the finite volume geometry vector for mpfa models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
#ifndef DUMUX_MIXEDDIMENSION_FACET_MPFA_GLOBALFVGEOMETRY_HH
#define DUMUX_MIXEDDIMENSION_FACET_MPFA_GLOBALFVGEOMETRY_HH

#include <dumux/discretization/cellcentered/mpfa/globalfvgeometry.hh>

namespace Dumux
{

template<class TypeTag>
class CCMpfaFacetGlobalFVGeometry : public CCMpfaGlobalFVGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>
{
    using ParentType = CCMpfaGlobalFVGeometry<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVGeometryCache)>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

public:
    //! Constructor
    CCMpfaFacetGlobalFVGeometry(const GridView gridView) : ParentType(gridView) {}

    //! After the update we have to initialize the coupling manager/mapper
    //! in case the global flux variables cache is enabled. This is because
    //! for the calculation of the transmissibilities we need
    void update(Problem& problem)
    {
        ParentType::update(problem);

        // maybe initialize coupling manager/mapper
        if (GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache))
            problem.couplingManager().init();
    }
};

} // end namespace

#endif
