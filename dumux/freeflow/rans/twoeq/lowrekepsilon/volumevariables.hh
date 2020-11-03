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
 * \ingroup LowReKEpsilonModel
 *
 * \copydoc Dumux::LowReKEpsilonVolumeVariables
 */
#ifndef DUMUX_LOWREKEPSILON_VOLUME_VARIABLES_HH
#define DUMUX_LOWREKEPSILON_VOLUME_VARIABLES_HH

#include <dumux/common/parameters.hh>
#include <dumux/freeflow/rans/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup LowReKEpsilonModel
 * \brief Volume variables for the isothermal single-phase low-Re k-epsilons model.
 */
template <class Traits, class NSVolumeVariables>
class LowReKEpsilonVolumeVariables
:  public RANSVolumeVariables<Traits, NSVolumeVariables>
{
    using RANSParentType = RANSVolumeVariables<Traits, NSVolumeVariables>;

    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    //! export the indices type
    using Indices = typename Traits::ModelTraits::Indices;

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
        RANSParentType::updateNavierStokesVolVars(elemSol, problem, element, scv);
        updateRANSProperties(elemSol, problem, element, scv);
    }

    /*!
     * \brief Update all turbulent quantities for a given control volume
     *
     * Wall and roughness related quantities are stored. Eddy viscosity is set.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void updateRANSProperties(const ElementSolution &elemSol,
                              const Problem &problem,
                              const Element &element,
                              const SubControlVolume& scv)
    {
        RANSParentType::updateRANSProperties(elemSol, problem, element, scv);
        turbulentKineticEnergy_ = elemSol[0][Indices::turbulentKineticEnergyIdx];
        dissipationTilde_ = elemSol[0][Indices::dissipationIdx];
        storedDissipationTilde_ = problem.storedDissipationTilde(RANSParentType::elementIdx());
        storedTurbulentKineticEnergy_ = problem.storedTurbulentKineticEnergy(RANSParentType::elementIdx());
        stressTensorScalarProduct_ = problem.stressTensorScalarProduct(RANSParentType::elementIdx());
        if (problem.useStoredEddyViscosity())
            RANSParentType::setDynamicEddyViscosity_(problem.storedDynamicEddyViscosity(RANSParentType::elementIdx()));
        else
            RANSParentType::setDynamicEddyViscosity_(calculateEddyViscosity());
        RANSParentType::calculateEddyDiffusivity(problem);
        RANSParentType::calculateEddyThermalConductivity(problem);
    }

    /*!
     * \brief Returns the dynamic eddy viscosity \f$\mathrm{[Pa s]}\f$.
     */
    Scalar calculateEddyViscosity()
    {
        return cMu() * fMu() * turbulentKineticEnergy() * turbulentKineticEnergy()
               / dissipationTilde() *  RANSParentType::density();
    }

    /*!
     * \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
     */
    Scalar turbulentKineticEnergy() const
    { return turbulentKineticEnergy_; }

    /*!
     * \brief Returns an effective dissipation \f$ m^2/s^3 \f$
     */
    Scalar dissipationTilde() const
    { return dissipationTilde_; }

    /*!
     * \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
     */
    Scalar storedTurbulentKineticEnergy() const
    { return storedTurbulentKineticEnergy_; }

    /*!
     * \brief Returns an effective dissipation \f$ m^2/s^3 \f$
     */
    Scalar storedDissipationTilde() const
    { return storedDissipationTilde_; }

    /*!
     * \brief Returns the scalar product of the stress tensor
     */
    Scalar stressTensorScalarProduct() const
    { return stressTensorScalarProduct_; }

    //! \brief Returns the \f$ Re_\textrm{T} \f$ value
    const Scalar reT() const
    {
        return turbulentKineticEnergy() * turbulentKineticEnergy()
               / RANSParentType::kinematicViscosity() / dissipationTilde();
    }

    //! \brief Returns the \f$ Re_\textrm{y} \f$ value
    const Scalar reY() const
    {
        using std::sqrt;
        return sqrt(turbulentKineticEnergy()) * RANSParentType::wallDistance()
               / RANSParentType::kinematicViscosity();
    }

    //! \brief Returns the \f$ C_\mu \f$ constant
    const Scalar cMu() const
    { return 0.09; }

    //! \brief Returns the \f$ \sigma_\textrm{k} \f$ constant
    const Scalar sigmaK() const
    { return 1.0; }

    //! \brief Returns the \f$ \sigma_\varepsilon \f$ constant
    const Scalar sigmaEpsilon() const
    { return 1.3; }

    //! \brief Returns the \f$ C_{1\tilde{\varepsilon}}  \f$ constant
    const Scalar cOneEpsilon() const
    { return 1.35; }

    //! \brief Returns the \f$ C_{2\tilde{\varepsilon}} \f$ constant
    const Scalar cTwoEpsilon() const
    { return 1.8; }

    //! \brief Returns the \f$ D \f$ value
    const Scalar dValue() const
    {
        return 2.0 * RANSParentType::kinematicViscosity() * turbulentKineticEnergy()
               / RANSParentType::wallDistance() / RANSParentType::wallDistance();
    }

    //! \brief Returns the \f$ f_\mu \f$ value
    const Scalar fMu() const
    {
        using std::exp;
        return 1.0 - exp(-0.0115 * RANSParentType::yPlus());
    }

    //! \brief Returns the \f$ f_1 \f$ value
    const Scalar fOne() const
    { return 1.0; }

    //! \brief Returns the \f$ f_2 \f$ value
    const Scalar fTwo() const
    {
        using std::exp;
        return 1.0 - 0.22 * exp(-1.0 * (reT() * reT() / 6.0 / 6.0));
    }

    //! \brief Returns the \f$ E \f$ value
    const Scalar eValue() const
    {
        using std::exp;
        return -2.0 * RANSParentType::kinematicViscosity() * dissipationTilde()
                / RANSParentType::wallDistance() / RANSParentType::wallDistance()
                * exp(-0.5 * RANSParentType::yPlus());
    }

protected:
    Scalar turbulentKineticEnergy_ = 0.0;
    Scalar dissipationTilde_ = 0.0;
    Scalar storedTurbulentKineticEnergy_ = 0.0;
    Scalar storedDissipationTilde_ = 0.0;
    Scalar stressTensorScalarProduct_ = 0.0;
};

} // end namespace Dumux

#endif
