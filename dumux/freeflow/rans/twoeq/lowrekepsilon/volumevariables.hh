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
 * \ingroup LowReKEpsilonModel
 *
 * \copydoc Dumux::LowReKEpsilonVolumeVariables
 */
#ifndef DUMUX_LOWREKEPSILON_VOLUME_VARIABLES_HH
#define DUMUX_LOWREKEPSILON_VOLUME_VARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/freeflow/rans/volumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup LowReKEpsilonModel
 * \brief Volume variables for the isothermal single-phase low-Re k-epsilons model.
 */
template <class Traits, class NSVolumeVariables>
class LowReKEpsilonVolumeVariables
:  public RANSVolumeVariables< Traits, LowReKEpsilonVolumeVariables<Traits, NSVolumeVariables> >
,  public NSVolumeVariables
{
    using ThisType = LowReKEpsilonVolumeVariables<Traits, NSVolumeVariables>;
    using RANSParentType = RANSVolumeVariables<Traits, ThisType>;
    using NavierStokesParentType = NSVolumeVariables;

    using Scalar = typename Traits::PrimaryVariables::value_type;

    static constexpr bool enableEnergyBalance = Traits::ModelTraits::enableEnergyBalance();
    static constexpr int fluidSystemPhaseIdx = Traits::ModelTraits::Indices::fluidSystemPhaseIdx;

public:
    //! export the underlying fluid system
    using FluidSystem = typename Traits::FluidSystem;
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
        NavierStokesParentType::update(elemSol, problem, element, scv);
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
        storedDissipationTilde_ = problem.storedDissipationTilde_[RANSParentType::elementID()];
        storedTurbulentKineticEnergy_ = problem.storedTurbulentKineticEnergy_[RANSParentType::elementID()];
        stressTensorScalarProduct_ = problem.stressTensorScalarProduct_[RANSParentType::elementID()];
        if (problem.useStoredEddyViscosity_)
            dynamicEddyViscosity_ = problem.storedDynamicEddyViscosity_[RANSParentType::elementID()];
        else
            dynamicEddyViscosity_ = calculateEddyViscosity();
        calculateEddyDiffusivity(problem);
    }

    /*!
     * \brief Return the dynamic eddy viscosity \f$\mathrm{[Pa s]}\f$ of the flow
     */
    Scalar dynamicEddyViscosity() const
    { return dynamicEddyViscosity_; }

    /*!
     * \brief Return the effective dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar effectiveViscosity() const
    { return NavierStokesParentType::viscosity() + dynamicEddyViscosity(); }

    /*!
     * \brief Returns the effective thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
     *        of the fluid-flow in the sub-control volume.
     */
    template<bool eB = enableEnergyBalance, typename std::enable_if_t<eB, int> = 0>
    Scalar effectiveThermalConductivity() const
    {
        return NavierStokesParentType::thermalConductivity()
               + RANSParentType::eddyThermalConductivity();
    }

    /*!
     * \brief Returns the dynamic eddy viscosity \f$\mathrm{[Pa s]}\f$.
     */
    Scalar calculateEddyViscosity()
    {
        return cMu() * fMu() * turbulentKineticEnergy() * turbulentKineticEnergy()
               / dissipationTilde() *  NavierStokesParentType::density();
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
        eddyDiffusivity_ = RANSParentType::kinematicEddyViscosity() / turbulentSchmidtNumber;
    }

    /*!
     * \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
     */
    Scalar turbulentKineticEnergy() const
    {
        return turbulentKineticEnergy_;
    }

    /*!
     * \brief Returns an effective dissipation \f$ m^2/s^3 \f$
     */
    Scalar dissipationTilde() const
    {
        return dissipationTilde_;
    }

    /*!
     * \brief Returns the turbulent kinetic energy \f$ m^2/s^2 \f$
     */
    Scalar storedTurbulentKineticEnergy() const
    {
        return storedTurbulentKineticEnergy_;
    }

    /*!
     * \brief Returns an effective dissipation \f$ m^2/s^3 \f$
     */
    Scalar storedDissipationTilde() const
    {
        return storedDissipationTilde_;
    }

    /*!
     * \brief Returns the scalar product of the stress tensor
     */
    Scalar stressTensorScalarProduct() const
    {
        return stressTensorScalarProduct_;
    }

    //! \brief Returns the \$f Re_\textrm{T} \$f value
    const Scalar reT() const
    {
        return turbulentKineticEnergy() * turbulentKineticEnergy()
               / RANSParentType::kinematicViscosity() / dissipationTilde();
    }

    //! \brief Returns the \$f Re_\textrm{y} \$f value
    const Scalar reY() const
    {
        using std::sqrt;
        return sqrt(turbulentKineticEnergy()) * RANSParentType::wallDistance()
               / RANSParentType::kinematicViscosity();
    }

    //! \brief Returns the \$f C_\mu \$f constant
    const Scalar cMu() const
    { return 0.09; }

    //! \brief Returns the \$f \sigma_\textrm{k} \$f constant
    const Scalar sigmaK() const
    { return 1.0; }

    //! \brief Returns the \$f \sigma_\varepsilon \$f constant
    const Scalar sigmaEpsilon() const
    { return 1.3; }

    //! \brief Returns the \$f C_{1\tilde{\varepsilon}}  \$f constant
    const Scalar cOneEpsilon() const
    { return 1.35; }

    //! \brief Returns the \$f C_{2\tilde{\varepsilon}} \$f constant
    const Scalar cTwoEpsilon() const
    { return 1.8; }

    //! \brief Returns the \$f D \$f value
    const Scalar dValue() const
    {
        return 2.0 * RANSParentType::kinematicViscosity() * turbulentKineticEnergy()
               / RANSParentType::wallDistance() / RANSParentType::wallDistance();
    }

    //! \brief Returns the \$f f_\mu \$f value
    const Scalar fMu() const
    {
        using std::exp;
        using std::pow;
        return 1.0 - exp(-0.0115 * RANSParentType::yPlus());
    }

    //! \brief Returns the \$f f_1 \$f value
    const Scalar fOne() const
    { return 1.0; }

    //! \brief Returns the \$f f_2 \$f value
    const Scalar fTwo() const
    {
        using std::exp;
        return 1.0 - 0.22 * exp(-1.0 * (reT() * reT() / 6.0 / 6.0));
    }

    //! \brief Returns the \$f E \$f value
    const Scalar eValue() const
    {
        using std::exp;
        return -2.0 * RANSParentType::kinematicViscosity() * dissipationTilde()
                / RANSParentType::wallDistance() / RANSParentType::wallDistance()
                * exp(-0.5 * RANSParentType::yPlus());
    }

    /*!
     * \brief Returns the eddy diffusivity \f$\mathrm{[m^2/s]}\f$
     */
    Scalar eddyDiffusivity() const
    { return eddyDiffusivity_; }

     /*!
     * \brief Returns the effective diffusion coefficient \f$\mathrm{[m^2/s]}\f$
     *
     * \param compIIdx the index of the component which diffusive
     * \param compJIdx the index of the component with respect to which compIIdx diffuses
     */
    Scalar effectiveDiffusivity(int compIIdx, int compJIdx = fluidSystemPhaseIdx) const
    {
        return NavierStokesParentType::diffusionCoefficient(compIIdx, compJIdx) + eddyDiffusivity();
    }

protected:
    Scalar dynamicEddyViscosity_;
    Scalar eddyDiffusivity_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipationTilde_;
    Scalar storedTurbulentKineticEnergy_;
    Scalar storedDissipationTilde_;
    Scalar stressTensorScalarProduct_;
};

}

#endif
