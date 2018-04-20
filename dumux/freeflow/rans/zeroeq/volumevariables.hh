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

#include <dune/common/exceptions.hh>
#include <dumux/freeflow/rans/volumevariables.hh>
#include "models.hh"

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
    using Indices = typename Traits::ModelTraits::Indices;

    static constexpr bool enableEnergyBalance = Traits::ModelTraits::enableEnergyBalance();



public:
    //! export the underlying fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! export the fluid state type
    using FluidState = typename Traits::FluidState;

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
        RANSParentType::updateRANSProperties(elemSol, problem, element, scv);
    }

    /*!
     * \brief Return the effective dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar effectiveViscosity() const
    {
        return NavierStokesParentType::viscosity()
               + RANSParentType::dynamicEddyViscosity();
    }

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
     */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    Scalar calculateEddyViscosity(const ElementSolution &elemSol,
                                  const Problem &problem,
                                  const Element &element,
                                  const SubControlVolume& scv)
    {
        additionalRoughnessLength_ = problem.additionalRoughnessLength_[RANSParentType::elementID()];
        yPlusRough_ = wallDistanceRough() * RANSParentType::uStar() / RANSParentType::kinematicViscosity();

        using std::abs;
        using std::exp;
        using std::sqrt;
        Scalar kinematicEddyViscosity = 0.0;
        unsigned int flowNormalAxis = problem.flowNormalAxis_[RANSParentType::elementID()];
        unsigned int wallNormalAxis = problem.wallNormalAxis_[RANSParentType::elementID()];
        Scalar velGrad = abs(RANSParentType::velocityGradients()[flowNormalAxis][wallNormalAxis]);

        if (problem.eddyViscosityModel_ == EddyViscosityModels::none)
        {
            // kinematicEddyViscosity = 0.0
        }
        else if (problem.eddyViscosityModel_ == EddyViscosityModels::prandtl)
        {
            Scalar mixingLength = problem.karmanConstant() * wallDistanceRough();
            kinematicEddyViscosity = mixingLength * mixingLength * velGrad;
        }
        else if (problem.eddyViscosityModel_ == EddyViscosityModels::modifiedVanDriest)
        {
            Scalar mixingLength = problem.karmanConstant() * wallDistanceRough()
                                  * (1.0 - exp(-yPlusRough() / 26.0))
                                  / sqrt(1.0 - exp(-0.26 * yPlusRough()));
            kinematicEddyViscosity = mixingLength * mixingLength * velGrad;
        }
        else if (problem.eddyViscosityModel_ == EddyViscosityModels::baldwinLomax)
        {
            kinematicEddyViscosity = problem.kinematicEddyViscosity_[RANSParentType::elementID()];
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented,
                       "This eddy viscosity model is not implemented: " << problem.eddyViscosityModel_);
        }

        return kinematicEddyViscosity * NavierStokesParentType::density();
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
     * \brief Return the fluid state of the control volume.
     */
    const FluidState& fluidState() const
    { return NavierStokesParentType::fluidState_; }

protected:

    Scalar additionalRoughnessLength_;
    Scalar yPlusRough_;
};

} // end namespace Dumux

#endif // DUMUX_ZEROEQ_VOLUME_VARIABLES_HH
