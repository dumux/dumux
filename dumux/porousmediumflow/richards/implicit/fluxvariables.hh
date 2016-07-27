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
 * \brief Overwrites the volumeFlux() function from the ImplicitDarcyFluxVariables.
 *
 */
#ifndef DUMUX_RICHARDS_FLUX_VARIABLES_HH
#define DUMUX_RICHARDS_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include "properties.hh"

namespace Dumux
{


/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitFluxVariables
 * \brief Overwrites the volumeFlux() function from the ImplicitDarcyFluxVariables.
 */
template <class TypeTag>
class RichardsFluxVariables : public ImplicitDarcyFluxVariables<TypeTag>
{
    friend class ImplicitDarcyFluxVariables<TypeTag>; // be friends with parent
    typedef ImplicitDarcyFluxVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { nPhaseIdx = Indices::nPhaseIdx} ;

public:
    /*!
     * \brief Return the volumetric flux over a face of a given phase.
     *
     *        This is the calculated velocity multiplied by the unit normal
     *        and the area of the face.
     *        face().normal
     *        has already the magnitude of the area.
     *        For the Richards model the velocity of the non-wetting phase
     *        is set to zero.
     *
     * \param phaseIdx index of the phase
     */
    Scalar volumeFlux(const unsigned int phaseIdx) const
    {
      if(phaseIdx == nPhaseIdx)
          return 0.;
      else
          return ParentType::volumeFlux(phaseIdx);
    }
};

} // end namespace Dumux

#endif
