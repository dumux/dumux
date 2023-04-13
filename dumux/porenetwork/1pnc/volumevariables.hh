// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
* \file
* \ingroup PNMOnePNCModel
* \copydoc Dumux::PoreNetwork::OnePNCVolumeVariables
*/

#ifndef DUMUX_PNM_1PNC_VOLUME_VARIABLES_HH
#define DUMUX_PNM_1PNC_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/1pnc/volumevariables.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMOnePNCModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the one-phase, n-component model.
 *
 * \note The default value for the phase index given in the fluid property interfaces is not used,
 *       but is only here to enable calling these functions without handing in a phase index
 *       (as in a single-phasic context there is only one phase).
 */
template <class Traits>
class OnePNCVolumeVariables
: public Dumux::OnePNCVolumeVariables<Traits>
{
    using ParentType = Dumux::OnePNCVolumeVariables<Traits>;
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

        poreInscribedRadius_ = problem.spatialParams().poreInscribedRadius(element, scv, elemSol);
        poreVolume_ = problem.gridGeometry().poreVolume(scv.dofIndex()) * this->porosity();
    }

    /*!
     * \brief Returns the pore's inscribed radius.
     */
    Scalar poreInscribedRadius() const
    { return poreInscribedRadius_; }

    /*!
     * \brief Returns the pore volume. // TODO should this be a fraction only?
     */
    Scalar poreVolume() const
    { return poreVolume_; }

protected:
    Scalar poreInscribedRadius_;
    Scalar poreVolume_;
};

} // end namespace Dumux::PoreNetwork

#endif
