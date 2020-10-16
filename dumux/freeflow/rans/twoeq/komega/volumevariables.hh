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
 * \ingroup KOmegaEqModel
 *
 * \copydoc Dumux::KOmegaVolumeVariables
 */
#ifndef DUMUX_KOMEGA_VOLUME_VARIABLES_HH
#define DUMUX_KOMEGA_VOLUME_VARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/freeflow/rans/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup KOmegaEqModel
 * \brief Volume variables for the isothermal single-phase k-omega 2-Eq model.
 */
template <class Traits, class NSVolumeVariables>
class KOmegaVolumeVariables
:  public RANSVolumeVariables<Traits, NSVolumeVariables>
{
    using RANSParentType = RANSVolumeVariables<Traits, NSVolumeVariables>;

    using Scalar = typename Traits::PrimaryVariables::value_type;
    using DimVector = Dune::FieldVector<Scalar, Traits::ModelTraits::dim()>;

public:
    //! export the indices type
    using Indices = typename Traits::ModelTraits::Indices;

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
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
        betaOmega_ = problem.betaOmega();
        turbulentKineticEnergy_ = elemSol[0][Indices::turbulentKineticEnergyIdx];
        dissipation_ = elemSol[0][Indices::dissipationIdx];
        storedDissipation_ = problem.storedDissipation(RANSParentType::elementIdx());
        storedTurbulentKineticEnergy_ = problem.storedTurbulentKineticEnergy(RANSParentType::elementIdx());
        storedDissipationGradient_ = problem.storedDissipationGradient(RANSParentType::elementIdx());
        storedTurbulentKineticEnergyGradient_ = problem.storedTurbulentKineticEnergyGradient(RANSParentType::elementIdx());
        stressTensorScalarProduct_ = problem.stressTensorScalarProduct(RANSParentType::elementIdx());
        if (problem.useStoredEddyViscosity())
            RANSParentType::setDynamicEddyViscosity_(problem.storedDynamicEddyViscosity(RANSParentType::elementIdx()));
        else
            RANSParentType::setDynamicEddyViscosity_(calculateEddyViscosity(problem));
        RANSParentType::calculateEddyDiffusivity(problem);
        RANSParentType::calculateEddyThermalConductivity(problem);
    }

    /*!
     * \brief Returns the dynamic eddy viscosity \f$\mathrm{[Pa s]}\f$.
     */
    template<class Problem>
    Scalar calculateEddyViscosity(const Problem& problem)
    {
        using std::sqrt;
        using std::max;

        // use the Dissipation limiter proposed in wilcox2008
        Scalar limitiedDissipation = std::numeric_limits<Scalar>::min();
        static const auto enableKOmegaDissipationLimiter
            = getParamFromGroup<bool>(problem.paramGroup(), "KOmega.EnableDissipationLimiter", true);
        if (enableKOmegaDissipationLimiter)
        {
            limitiedDissipation = (7.0 / 8.0) * sqrt(2.0 * stressTensorScalarProduct() / betaK());
        }
        return turbulentKineticEnergy() / max(dissipation(), limitiedDissipation)
               * RANSParentType::density();
    }

    //! \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
    Scalar turbulentKineticEnergy() const
    { return turbulentKineticEnergy_; }

    //! \brief Returns an effective dissipation \f$ m^2/s^3 \f$
    Scalar dissipation() const
    { return dissipation_; }

    //! \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
    Scalar storedTurbulentKineticEnergy() const
    { return storedTurbulentKineticEnergy_; }

    //! \brief Returns an effective dissipation \f$ m^2/s^3 \f$
    Scalar storedDissipation() const
    { return storedDissipation_; }

    //! \brief Returns the gradient of the turbulent kinetic energy \f$ m^2/s^2 \f$
    DimVector storedTurbulentKineticEnergyGradient() const
    { return storedTurbulentKineticEnergyGradient_; }

    //! \brief Returns the gradient of the effective dissipation \f$ m^2/s^3 \f$
    DimVector storedDissipationGradient() const
    { return storedDissipationGradient_; }

    //! \brief Returns the scalar product of the stress tensor
    Scalar stressTensorScalarProduct() const
    { return stressTensorScalarProduct_; }

    //! \brief Returns the \f$ \alpha \f$ value
    const Scalar alpha() const
    { return 0.52; }

    //! \brief Returns the \f$ \sigma_k \f$ constant
    const Scalar sigmaK() const
    { return 0.6; }

    //! \brief Returns the \f$ \sigma_{\omega} \f$ constant
    const Scalar sigmaOmega() const
    { return 0.5; }

    //! \brief Returns the \f$ \beta_k \f$ constant
    const Scalar betaK() const
    { return 0.09; }

    //! \brief Returns the \f$ \beta_\omega \f$ constant
    const Scalar betaOmega() const
    { return betaOmega_; }

protected:
    Scalar betaOmega_ = 0.0;
    Scalar dissipation_ = 0.0;
    Scalar turbulentKineticEnergy_ = 0.0;
    Scalar storedDissipation_ = 0.0;
    DimVector storedDissipationGradient_ = DimVector(0.0);
    Scalar storedTurbulentKineticEnergy_ = 0.0;
    DimVector storedTurbulentKineticEnergyGradient_ = DimVector(0.0);
    Scalar stressTensorScalarProduct_ = 0.0;
};

} // end namespace Dumux

#endif
