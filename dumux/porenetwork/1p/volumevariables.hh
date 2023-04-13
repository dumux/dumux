// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
* \file
* \ingroup PNMOnePModel
* \copydoc Dumux::PoreNetwork::OnePVolumeVariables
*/

#ifndef DUMUX_PNM_1P_VOLUME_VARIABLES_HH
#define DUMUX_PNM_1P_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/1p/volumevariables.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMOnePModel
 * \brief Contains the quantities which are constant within a
 *        finite volume (the pore body) in the one-phase model.
 *
 * \tparam Traits Class encapsulating types to be used by the volVars
 */
template<class Traits>
class OnePVolumeVariables : public Dumux::OnePVolumeVariables<Traits>
{
    using ParentType = Dumux::OnePVolumeVariables<Traits>;
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
        inscribedPoreRadius_ = problem.spatialParams().poreInscribedRadius(element, scv, elemSol);
        poreVolume_ = problem.gridGeometry().poreVolume(scv.dofIndex()) * this->porosity();
    }

    /*!
     * \brief Returns the pore's inscribed radius.
     */
    Scalar poreInscribedRadius() const
    { return inscribedPoreRadius_; }

    /*!
     * \brief Returns the pore volume. // TODO should this be a fraction only?
     */
    Scalar poreVolume() const
    { return poreVolume_; }

protected:
    Scalar inscribedPoreRadius_;
    Scalar poreVolume_;
};

} // end namespace Dumux::PoreNetwork

#endif
