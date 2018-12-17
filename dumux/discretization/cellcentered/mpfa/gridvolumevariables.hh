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
 * \ingroup CCMpfaDiscretization
 * \brief The grid volume variables class for cell centered mpfa models
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_GRID_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_GRID_VOLUMEVARIABLES_HH

#include <dumux/discretization/cellcentered/mpfa/elementvolumevariables.hh>
#include <dumux/discretization/cellcentered/gridvolumevariables.hh>

namespace Dumux {

template<class P, class VV>
struct CCMpfaDefaultGridVolumeVariablesTraits
{
    using Problem = P;
    using VolumeVariables = VV;

    template<class GridVolumeVariables, bool cachingEnabled>
    using LocalView = CCMpfaElementVolumeVariables<GridVolumeVariables, cachingEnabled>;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Base class for the grid volume variables
 * \note This class has a cached version and a non-cached version
 * \tparam Problem the type of problem we are solving
 * \tparam VolumeVariables the type of volume variables we are using for the model
 * \tparam Traits the traits class injecting the problem, volVar and elemVolVars type
 * \tparam cachingEnabled if the cache is enabled
 */
template<class Problem,
         class VolumeVariables,
         bool cachingEnabled = false,
         class Traits = CCMpfaDefaultGridVolumeVariablesTraits<Problem, VolumeVariables> >
class CCMpfaGridVolumeVariables : public CCGridVolumeVariables<Traits, cachingEnabled>
{
public:
    using ParentType = CCGridVolumeVariables<Traits, cachingEnabled>;
    using ParentType::ParentType;
};

} // end namespace Dumux

#endif
