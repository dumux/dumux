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
 * \ingroup KEpsilonModel
 * \copydoc Dumux::KEpsilonVolumeVariables
 */
#ifndef DUMUX_KEPSILON_VOLUME_VARIABLES_HH
#define DUMUX_KEPSILON_VOLUME_VARIABLES_HH

#include <dumux/common/parameters.hh>
#include <dumux/freeflow/rans/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup KEpsilonModel
 * \brief Volume variables for the isothermal single-phase k-epsilon model.
 */
template <class Traits, class NSVolumeVariables>
class KEpsilonVolumeVariables
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
        isMatchingPoint_ = problem.isMatchingPoint(RANSParentType::elementIdx());
        inNearWallRegion_ = problem.inNearWallRegion(RANSParentType::elementIdx());
        turbulentKineticEnergy_ = elemSol[0][Indices::turbulentKineticEnergyIdx];
        dissipation_ = elemSol[0][Indices::dissipationIdx];
        storedDissipation_ = problem.storedDissipation(RANSParentType::elementIdx());
        storedTurbulentKineticEnergy_ = problem.storedTurbulentKineticEnergy(RANSParentType::elementIdx());
        stressTensorScalarProduct_ = problem.stressTensorScalarProduct(RANSParentType::elementIdx());
        const Scalar uStarNominal = problem.uStarNominal(RANSParentType::elementIdx());
        const auto flowDirectionAxis = problem.flowDirectionAxis(RANSParentType::elementIdx());
        yPlusNominal_ = RANSParentType::wallDistance() * uStarNominal / problem.kinematicViscosity(RANSParentType::elementIdx());
        uPlusNominal_ = RANSParentType::ccVelocityVector()[flowDirectionAxis] / uStarNominal;
        cMu_ = problem.cMu();
        if (problem.useStoredEddyViscosity())
            RANSParentType::setDynamicEddyViscosity_(problem.storedDynamicEddyViscosity(RANSParentType::elementIdx()));
        else
            RANSParentType::setDynamicEddyViscosity_(calculateEddyViscosity());
        if (inNearWallRegion_)
        {
            RANSParentType::setDynamicEddyViscosity_(problem.zeroEqDynamicEddyViscosity(RANSParentType::elementIdx()));
        }
        RANSParentType::calculateEddyDiffusivity(problem);
        RANSParentType::calculateEddyThermalConductivity(problem);
    }

    /*!
     * \brief Returns the dynamic eddy viscosity \f$\mathrm{[Pa s]}\f$.
     */
    Scalar calculateEddyViscosity()
    {
        return cMu() * turbulentKineticEnergy() * turbulentKineticEnergy()
               / dissipation() *  RANSParentType::density();
    }

    /*!
     * \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
     */
    Scalar turbulentKineticEnergy() const
    { return turbulentKineticEnergy_; }

    /*!
     * \brief Returns an effective dissipation \f$ m^2/s^3 \f$
     */
    Scalar dissipation() const
    { return dissipation_; }

    /*!
     * \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
     */
    Scalar storedTurbulentKineticEnergy() const
    { return storedTurbulentKineticEnergy_; }

    /*!
     * \brief Returns an effective dissipation \f$ m^2/s^3 \f$
     */
    Scalar storedDissipation() const
    { return storedDissipation_; }

    /*!
     * \brief Returns the scalar product of the stress tensor
     */
    Scalar stressTensorScalarProduct() const
    { return stressTensorScalarProduct_; }

    /*
     * \brief Returns if an element is located in the near-wall region
     */
    bool inNearWallRegion() const
    { return inNearWallRegion_; }

    /*!
     * \brief Returns if an element is the matching point
     */
    bool isMatchingPoint() const
    { return isMatchingPoint_; }

    //! \brief Returns the \f$ C_{\mu} \f$ constant
    const Scalar cMu() const
    { return cMu_; }

    //! \brief Returns the \f$ \sigma_{\textrm{k}} \f$ constant
    const Scalar sigmaK() const
    { return 1.0; }

    //! \brief Returns the \f$ \sigma_{\varepsilon} \f$ constant
    const Scalar sigmaEpsilon() const
    { return 1.3; }

    //! \brief Returns the \f$ C_{1\varepsilon}  \f$ constant
    const Scalar cOneEpsilon() const
    { return 1.44; }

    //! \brief Returns the \f$ C_{2\varepsilon} \f$ constant
    const Scalar cTwoEpsilon() const
    { return 1.92; }

    //! \brief Returns the nominal dimensionless wall distance \f$\mathrm{[-]}\f$.
    Scalar yPlusNominal() const
    { return yPlusNominal_; }

    //! \brief Return the nominal dimensionless velocity \f$\mathrm{[-]}\f$.
    Scalar uPlusNominal() const
    { return uPlusNominal_; }

protected:
    Scalar turbulentKineticEnergy_ = 0.0;
    Scalar dissipation_ = 0.0;
    Scalar storedTurbulentKineticEnergy_ = 0.0;
    Scalar storedDissipation_ = 0.0;
    Scalar stressTensorScalarProduct_ = 0.0;
    Scalar yPlusNominal_ = 0.0;
    Scalar uPlusNominal_ = 0.0;
    Scalar cMu_ = 0.0;
    bool inNearWallRegion_ = false;
    bool isMatchingPoint_ = false;
};

} // end namespace Dumux

#endif
