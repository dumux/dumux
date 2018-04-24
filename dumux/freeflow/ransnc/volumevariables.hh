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
 * \ingroup RANSNCModel
 *
 * \copydoc Dumux::RANSNCVolumeVariables
 */
#ifndef DUMUX_RANS_NC_VOLUMEVARIABLES_HH
#define DUMUX_RANS_NC_VOLUMEVARIABLES_HH

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup RANSNCModel
 * \brief Volume variables for the single-phase, multi-component Reynolds-averaged Navier-Stokes model.
 */
template <class Traits, class RANSVolVars>
class RANSNCVolumeVariables : public RANSVolVars
{
    using ParentType = RANSVolVars;
    using Scalar = typename Traits::PrimaryVariables::value_type;

    static constexpr int fluidSystemPhaseIdx = Traits::ModelTraits::Indices::fluidSystemPhaseIdx;

public:

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);
        calculateEddyDiffusivity(problem);
    }

    /*!
     * \brief Calculates the eddy diffusivity \f$\mathrm{[m^2/s]}\f$ based
     *        on the kinematic eddy viscosity and the turbulent schmidt number
     */
    template<class Problem>
    void calculateEddyDiffusivity(const Problem& problem)
    {
        static const auto turbulentSchmidtNumber
            = getParamFromGroup<Scalar>(problem.paramGroup(),
                                        "RANS.TurbulentSchmidtNumber", 1.0);
        eddyDiffusivity_ = ParentType::kinematicEddyViscosity()
                           / turbulentSchmidtNumber;
    }

    /*!
     * \brief Returns the eddy diffusivity \f$\mathrm{[m^2/s]}\f$
     */
    Scalar eddyDiffusivity() const
    {
        return eddyDiffusivity_;
    }

    /*!
     * \brief Returns the effective diffusion coefficient \f$\mathrm{[m^2/s]}\f$
     *
     * \param compIIdx the index of the component which diffusive
     * \param compJIdx the index of the component with respect to which compIIdx diffuses
     */
    Scalar effectiveDiffusivity(int compIIdx, int compJIdx = fluidSystemPhaseIdx) const
    {
        return ParentType::diffusionCoefficient(compIIdx, fluidSystemPhaseIdx)
               + eddyDiffusivity();
    }

private:
    Scalar eddyDiffusivity_;
};

} // end namespace Dumux

#endif
