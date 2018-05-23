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
 * \ingroup KEpsilonModel
 *
 * \copydoc Dumux::KEpsilonVolumeVariables
 */
#ifndef DUMUX_KEPSILON_VOLUME_VARIABLES_HH
#define DUMUX_KEPSILON_VOLUME_VARIABLES_HH

#include <dumux/common/parameters.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/freeflow/rans/volumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup KEpsilonModel
 * \brief Volume variables for the isothermal single-phase k-epsilon model.
 */
template <class Traits, class NSVolumeVariables>
class KEpsilonVolumeVariables
:  public RANSVolumeVariables< Traits, KEpsilonVolumeVariables<Traits, NSVolumeVariables> >
,  public NSVolumeVariables
{
    using ThisType = KEpsilonVolumeVariables<Traits, NSVolumeVariables>;
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
        isMatchingPoint_ = problem.isMatchingPoint(RANSParentType::elementID());
        inNearWallRegion_ = problem.inNearWallRegion(RANSParentType::elementID());
        turbulentKineticEnergy_ = elemSol[0][Indices::turbulentKineticEnergyIdx];
        dissipation_ = elemSol[0][Indices::dissipationIdx];
        storedDissipation_ = problem.storedDissipation_[RANSParentType::elementID()];
        storedTurbulentKineticEnergy_ = problem.storedTurbulentKineticEnergy_[RANSParentType::elementID()];
        stressTensorScalarProduct_ = problem.stressTensorScalarProduct_[RANSParentType::elementID()];
        const Scalar uStarNominal = problem.uStarNominal(RANSParentType::elementID());
        const auto flowNormalAxis = problem.flowNormalAxis_[RANSParentType::elementID()];
        yPlusNominal_ = RANSParentType::wallDistance() * uStarNominal / problem.kinematicViscosity_[RANSParentType::elementID()];
        uPlusNominal_ = RANSParentType::velocity()[flowNormalAxis] / uStarNominal;
        if (problem.useStoredEddyViscosity_)
            dynamicEddyViscosity_ = problem.storedDynamicEddyViscosity_[RANSParentType::elementID()];
        else
            dynamicEddyViscosity_ = calculateEddyViscosity();
        if (inNearWallRegion_ && !isMatchingPoint_)
        {
            dynamicEddyViscosity_ = problem.zeroEqDynamicEddyViscosity_[RANSParentType::elementID()];
        }
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
        return cMu() * turbulentKineticEnergy() * turbulentKineticEnergy()
               / dissipation() *  NavierStokesParentType::density();
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
    Scalar dissipation() const
    {
        return dissipation_;
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
    Scalar storedDissipation() const
    {
        return storedDissipation_;
    }

    /*!
     * \brief Returns the scalar product of the stress tensor
     */
    Scalar stressTensorScalarProduct() const
    {
        return stressTensorScalarProduct_;
    }

    /*
     * \brief Returns if an element is located in the near-wall region
     */
    bool inNearWallRegion() const
    {
        return inNearWallRegion_;
    }

    /*!
     * \brief Returns if an element is the matching point
     */
    Scalar isMatchingPoint() const
    {
        return isMatchingPoint_;
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

    //! \brief Returns the \$f C_{1\varepsilon}  \$f constant
    const Scalar cOneEpsilon() const
    { return 1.44; }

    //! \brief Returns the \$f C_{2\varepsilon} \$f constant
    const Scalar cTwoEpsilon() const
    { return 1.92; }

    /*!
     * \brief Return the nominal dimensionless wall distance \f$\mathrm{[-]}\f$.
     */
    Scalar yPlusNominal() const
    { return yPlusNominal_; }

    /*!
     * \brief Return the nominal dimensionless velocity \f$\mathrm{[-]}\f$.
     */
    Scalar uPlusNominal() const
    { return uPlusNominal_; }

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
    Scalar dissipation_;
    Scalar storedTurbulentKineticEnergy_;
    Scalar storedDissipation_;
    Scalar stressTensorScalarProduct_;
    Scalar yPlusNominal_;
    Scalar uPlusNominal_;
    bool inNearWallRegion_;
    bool isMatchingPoint_;
};

}

#endif
