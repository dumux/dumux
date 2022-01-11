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
 * \ingroup TwoEqModel
 *
 * \copydoc Dumux::TwoEqVolumeVariables
 */
#ifndef DUMUX_TWOEQ_VOLUME_VARIABLES_HH
#define DUMUX_TWOEQ_VOLUME_VARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dune/common/fvector.hh>
#include <dumux/freeflow/navierstokes/turbulence/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup TwoEqModel
 * \brief Volume variables for the isothermal single-phase 2-Eq model.
 */
template <class Traits, class NSVolumeVariables, int dim>
class TwoEqVolumeVariables
:  public RANSVolumeVariables<Traits, NSVolumeVariables, dim>
{
    using RANSParentType = RANSVolumeVariables<Traits, NSVolumeVariables, dim>;

    using Scalar = typename Traits::PrimaryVariables::value_type;
    using DimVector = Dune::FieldVector<Scalar, dim>;

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
        RANSParentType::update(elemSol, problem, element, scv);
        updateRANSVariables(elemSol, problem, element, scv);
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
    void updateRANSVariables(const ElementSolution &elemSol,
                             const Problem &problem,
                             const Element &element,
                             const SubControlVolume& scv)
    {
        // Update the base Rans Volvars
        RANSParentType::updateRANSVariables(elemSol, problem, element, scv);

        // Save the primary variables
        turbulentKineticEnergy_ = elemSol[0][Indices::turbulentKineticEnergyIdx];
        dissipation_ = elemSol[0][Indices::dissipationIdx];
        twoEqTurbulenceModelName_ = problem.twoEqTurbulenceModelName();

        // Solve the eddy viscosity
        updateEddyViscositites();

        // Update the compositional and non-isothermal diffusive coeffs
        RANSParentType::calculateEddyDiffusivity();
        RANSParentType::calculateEddyThermalConductivity();
    }

    void updateEddyViscositites()
    {
        if (twoEqTurbulenceModelName() == "Wilcox")
        {
            using std::sqrt;
            using std::max;
            // use the Dissipation limiter proposed in wilcox2008
            Scalar limitiedDissipation = std::numeric_limits<Scalar>::min();
            static const auto enableKOmegaDissipationLimiter = getParam<bool>("RANS.EnableKOmegaDissipationLimiter", true);
            if (enableKOmegaDissipationLimiter)
                limitiedDissipation = (7.0 / 8.0) * sqrt(2.0 * this->stressTensorScalarProduct() / betaK());

            dynamicEddyViscosity_ = turbulentKineticEnergy() / max(dissipation(), limitiedDissipation) * RANSParentType::density();
        }
        else if (twoEqTurbulenceModelName() == "SST")
        {
            //TODO:
        }
        else if (twoEqTurbulenceModelName() == "BSL")
        {
            //TODO:
        }
        else if (twoEqTurbulenceModelName() == "LowReKEpsilon")
        {
            //TODO:
        }
        else if (twoEqTurbulenceModelName() == "KEpsilon")
        {
            //TODO:
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented, "This turbulence model " << twoEqTurbulenceModelName() <<
                                             " is not available, try a different two-eq model.");
        }
    }

    //! \brief Returns the dynamic Eddy Viscosity \f$\mathrm{[Pa s]}\f$.
    Scalar dynamicEddyViscosity() const
    { return dynamicEddyViscosity_; }

    //! \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
    Scalar turbulentKineticEnergy() const
    { return turbulentKineticEnergy_; }

    //! \brief Returns an effective dissipation \f$ m^2/s^3 \f$
    Scalar dissipation() const
    { return dissipation_; }

    Scalar dissipationOmega() const
    {
        if (twoEqTurbulenceModelName() == "SST" ||
            twoEqTurbulenceModelName() == "BSL" ||
            twoEqTurbulenceModelName() == "Wilcox" )
            return dissipation();
        else if (twoEqTurbulenceModelName() == "LowReKEpsilon" ||
                 twoEqTurbulenceModelName() == "KEpsilon")
            return dissipation() / turbulentKineticEnergy();
        else
            DUNE_THROW(Dune::NotImplemented, "This turbulence model " << twoEqTurbulenceModelName() <<
                                             " is not available, try a different two-eq model.");
    }

    Scalar dissipationEpsilon() const
    {
        if  (twoEqTurbulenceModelName() == "LowReKEpsilon" ||
             twoEqTurbulenceModelName() == "KEpsilon")
            return dissipation();
        else if(twoEqTurbulenceModelName() == "SST" ||
                twoEqTurbulenceModelName() == "BSL" ||
                twoEqTurbulenceModelName() == "Wilcox" )
            return dissipation() * turbulentKineticEnergy();
        else
            DUNE_THROW(Dune::NotImplemented, "This turbulence model " << twoEqTurbulenceModelName() <<
                                             " is not available, try a different two-eq model.");
    }
    //! \brief Returns the gradient of the turbulent kinetic energy \f$ m^2/s^2 \f$
    DimVector storedTurbulentKineticEnergyGradient() const
    { return DimVector(0.0); }// TODO: Figure this out

    //! \brief Returns the gradient of the effective dissipation \f$ m^2/s^3 \f$
    DimVector storedDissipationGradient() const
    { return DimVector(0.0); }// TODO: Figure this out

    const std::string twoEqTurbulenceModelName() const
    { return twoEqTurbulenceModelName_; }

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
    { return 0.0708; }

protected:
    Scalar dissipation_ = 0.0;
    Scalar turbulentKineticEnergy_ = 0.0;
    Scalar dynamicEddyViscosity_;
    std::string twoEqTurbulenceModelName_;
};

} // end namespace Dumux

#endif
