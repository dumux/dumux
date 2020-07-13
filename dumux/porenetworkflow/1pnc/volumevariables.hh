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
 * \ingroup OnePNCModel
 * \brief Quantities required by the single-phase, n-component box
 *        model defined on a vertex.
 */

#ifndef DUMUX_PNM_1PNC_VOLUME_VARIABLES_HH
#define DUMUX_PNM_1PNC_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/1pnc/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the one-phase, n-component model.
 *
 * \note The default value for the phase index given in the fluid property interfaces is not used,
 *       but is only here to enable calling these functions without handing in a phase index
 *       (as in a single-phasic context there is only one phase).
 */
template <class Traits>
class PNMOnePNCVolumeVariables
: public OnePNCVolumeVariables<Traits>
{
    using ParentType = OnePNCVolumeVariables<Traits>;
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:

    /*!
     * \brief Updates all quantities for a given control volume.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        poreRadius_ = problem.spatialParams().poreRadius(element, scv, elemSol);
        poreVolume_ = problem.gridGeometry().poreVolume(scv.dofIndex()) * this->porosity();
    }

    Scalar poreRadius() const
    { return poreRadius_; }

    Scalar poreVolume() const
    { return poreVolume_; }

protected:

    Scalar poreRadius_;
    Scalar poreVolume_;
};

} // end namespace Dumux

#endif
