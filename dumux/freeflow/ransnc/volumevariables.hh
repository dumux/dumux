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

#include <dumux/common/properties.hh>

#include <dumux/freeflow/rans/volumevariables.hh>
#include <dumux/freeflow/navierstokesnc/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup RANSNCModel
 * \brief Volume variables for the single-phase, multi-component Reynolds-averaged Navier-Stokes model.
 */
template <class TypeTag>
class RANSNCVolumeVariables
    : virtual public GET_PROP_TYPE(TypeTag, SinglePhaseVolumeVariables),
      virtual public NavierStokesNCVolumeVariables<TypeTag>
{
    using ParentTypeSinglePhase = typename GET_PROP_TYPE(TypeTag, SinglePhaseVolumeVariables);
    using ParentTypeCompositional = NavierStokesNCVolumeVariables<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr auto phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);

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
    template<class ElementSolution>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentTypeCompositional::update(elemSol, problem, element, scv);
        ParentTypeSinglePhase::updateRANSProperties(elemSol, problem, element, scv);
        calculateEddyDiffusivity();
    }

    /*!
     * \brief Calculates the eddy diffusivity \f$\mathrm{[m^2/s]}\f$ based
     *        on the kinematic eddy viscosity and the turbulent schmidt number
     */
    void calculateEddyDiffusivity()
    {
        static const auto turbulentSchmidtNumber
            = getParamFromGroup<Scalar>(GET_PROP_VALUE(TypeTag, ModelParameterGroup),
                                        "RANS.TurbulentSchmidtNumber", 1.0);
        eddyDiffusivity_ = ParentTypeSinglePhase::kinematicEddyViscosity()
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
     * \param compIdx the index of the component
     */
    Scalar effectiveDiffusivity(int compIdx) const
    {
        return ParentTypeCompositional::diffusionCoefficient(compIdx)
               + eddyDiffusivity();
    }

protected:
    Scalar eddyDiffusivity_;
};

} // end namespace Dumux

#endif
