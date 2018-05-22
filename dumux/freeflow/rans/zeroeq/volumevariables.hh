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
 * \ingroup ZeroEqModel
 *
 * \copydoc Dumux::ZeroEqVolumeVariables
 */
#ifndef DUMUX_ZEROEQ_VOLUME_VARIABLES_HH
#define DUMUX_ZEROEQ_VOLUME_VARIABLES_HH

#include <string>

#include <dune/common/exceptions.hh>
#include <dumux/freeflow/rans/volumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup ZeroEqModel
 * \brief Volume variables for the single-phase 0-Eq. model.
 */
template <class Traits, class NSVolumeVariables>
class ZeroEqVolumeVariables
:  public RANSVolumeVariables< Traits, ZeroEqVolumeVariables<Traits, NSVolumeVariables> >
,  public NSVolumeVariables
{
    using ThisType = ZeroEqVolumeVariables<Traits, NSVolumeVariables>;
    using RANSParentType = RANSVolumeVariables<Traits, ThisType>;
    using NavierStokesParentType = NSVolumeVariables;

    using Scalar = typename Traits::PrimaryVariables::value_type;

    static constexpr bool enableEnergyBalance = Traits::ModelTraits::enableEnergyBalance();
    static constexpr int fluidSystemPhaseIdx = Traits::ModelTraits::Indices::fluidSystemPhaseIdx;

public:
    //! export the underlying fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! export the fluid state type
    using FluidState = typename Traits::FluidState;
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
        additionalRoughnessLength_ = problem.additionalRoughnessLength_[RANSParentType::elementID()];
        yPlusRough_ = wallDistanceRough() * RANSParentType::uStar() / RANSParentType::kinematicViscosity();
        dynamicEddyViscosity_ = calculateEddyViscosity(elemSol, problem, element, scv, problem.eddyViscosityModel_);
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
     * \brief Calculate and set the dynamic eddy viscosity.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     * \param modelName The name of the used model
     */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    Scalar calculateEddyViscosity(const ElementSolution &elemSol,
                                  const Problem &problem,
                                  const Element &element,
                                  const SubControlVolume& scv,
                                  const std::string modelName)
    {
        using std::abs;
        using std::exp;
        using std::sqrt;
        Scalar kinematicEddyViscosity = 0.0;
        unsigned int flowNormalAxis = problem.flowNormalAxis_[RANSParentType::elementID()];
        unsigned int wallNormalAxis = problem.wallNormalAxis_[RANSParentType::elementID()];
        Scalar velGrad = abs(RANSParentType::velocityGradients()[flowNormalAxis][wallNormalAxis]);

        if (modelName.compare("none") == 0)
        {
            // kinematicEddyViscosity = 0.0
        }
        else if (modelName.compare("prandtl") == 0)
        {
            Scalar mixingLength = problem.karmanConstant() * wallDistanceRough();
            kinematicEddyViscosity = mixingLength * mixingLength * velGrad;
        }
        else if (modelName.compare("vanDriest") == 0)
        {
            Scalar mixingLength = problem.karmanConstant() * wallDistanceRough()
                                  * (1.0 - exp(-yPlusRough() / 26.0))
                                  / sqrt(1.0 - exp(-0.26 * yPlusRough()));
            kinematicEddyViscosity = mixingLength * mixingLength * velGrad;
        }
        else if (modelName.compare("baldwinLomax") == 0)
        {
            kinematicEddyViscosity = problem.kinematicEddyViscosity_[RANSParentType::elementID()];
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented,
                       "The eddy viscosity model \"" << modelName << "\" is not implemented.");
        }

        return kinematicEddyViscosity * NavierStokesParentType::density();
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
     * \brief Return the wall distance \f$\mathrm{[m]}\f$ including an additional roughness length
     */
    Scalar wallDistanceRough() const
    { return RANSParentType::wallDistance() + additionalRoughnessLength_; }

    /*!
     * \brief Return the dimensionless wall distance \f$\mathrm{[-]}\f$  including an additional roughness length
     */
    Scalar yPlusRough() const
    { return yPlusRough_; }

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

    /*!
     * \brief Return the fluid state of the control volume.
     */
    const FluidState& fluidState() const
    { return NavierStokesParentType::fluidState_; }

protected:
    Scalar dynamicEddyViscosity_;
    Scalar eddyDiffusivity_;
    Scalar additionalRoughnessLength_;
    Scalar yPlusRough_;
};

} // end namespace Dumux

#endif // DUMUX_ZEROEQ_VOLUME_VARIABLES_HH
