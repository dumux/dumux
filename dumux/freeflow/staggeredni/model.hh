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
 *
 * \brief Base class for all models which use the one-phase,
 *        fully implicit model.
 *        Adaption of the fully implicit scheme to the one-phase flow model.
 */

#ifndef DUMUX_STAGGERED_NI_MODEL_HH
#define DUMUX_STAGGERED_NI_MODEL_HH

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup NavierStokesModel
 * \brief A single-phase, non-isothermal flow model using the fully implicit scheme.
 *
 * All equations are discretized using a staggered grid as spatial
 * and the implicit Euler method as time discretization.
 * The model supports compressible as well as incompressible fluids.
 */
template<class TypeTag >
class NavierStokesNonIsothermalModel : public GET_PROP_TYPE(TypeTag, IsothermalModel)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, IsothermalModel);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

public:

    void init(Problem& problem)
    {
        ParentType::init(problem);

        // add temperature to output
        auto& vtkOutputModule = problem.vtkOutputModule();
        vtkOutputModule.addPrimaryVariable("temperature", Indices::temperatureIdx);
    }
};
}

#include "propertydefaults.hh"

#endif
