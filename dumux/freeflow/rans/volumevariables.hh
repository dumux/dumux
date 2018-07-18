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
 * \ingroup RANSModel
 *
 * \copydoc Dumux::RANSVolumeVariables
 */
#ifndef DUMUX_RANS_VOLUME_VARIABLES_HH
#define DUMUX_RANS_VOLUME_VARIABLES_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/parameters.hh>

namespace Dumux
{
/*!
 * \ingroup RANSModel
 * \brief Volume variables for the isothermal single-phase Reynolds-Averaged Navier-Stokes models.
 */
template <class Traits, class NSVolumeVariables>
class RANSVolumeVariables
: public NSVolumeVariables
{
    using NavierStokesParentType = NSVolumeVariables;

    using Scalar = typename Traits::PrimaryVariables::value_type;

    enum { dimWorld = Traits::ModelTraits::dim() };
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    static constexpr bool enableEnergyBalance = Traits::ModelTraits::enableEnergyBalance();

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
    void updateNavierStokesVolVars(const ElementSolution &elemSol,
                                   const Problem &problem,
                                   const Element &element,
                                   const SubControlVolume& scv)
    {
        NavierStokesParentType::update(elemSol, problem, element, scv);
    }

    /*!
     * \brief Update all turbulent quantities for a given control volume
     *
     * Wall related quantities are stored and the calculateEddyViscosity(...)
     * function of the turbulence model implementation is called.
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
        using std::abs;
        using std::max;
        using std::sqrt;

        // calculate characteristic properties of the turbulent flow
        elementID_ = problem.fvGridGeometry().elementMapper().index(element);
        wallElementID_ = problem.wallElementID_[elementID_];
        wallDistance_ = problem.wallDistance_[elementID_];
        velocity_ = problem.velocity_[elementID_];
        velocityMaximum_ = problem.velocityMaximum_[wallElementID_];
        velocityGradients_ = problem.velocityGradients_[elementID_];
        const auto flowNormalAxis = problem.flowNormalAxis_[elementID_];
        const auto wallNormalAxis = problem.wallNormalAxis_[elementID_];
        uStar_ = sqrt(problem.kinematicViscosity_[wallElementID_]
                      * abs(problem.velocityGradients_[wallElementID_][flowNormalAxis][wallNormalAxis]));
        uStar_ = max(uStar_, 1e-10); // zero values lead to numerical problems in some turbulence models
        yPlus_ = wallDistance_ * uStar_ / problem.kinematicViscosity_[elementID_];
        uPlus_ = velocity_[flowNormalAxis] / uStar_;
        karmanConstant_ = problem.karmanConstant();
    }

    /*!
     * \brief Return the element ID of the control volume.
     */
    unsigned int elementID() const
    { return elementID_; }

    /*!
     * \brief Return the velocity vector \f$\mathrm{[m/s]}\f$ at the control volume center.
     */
    DimVector velocity() const
    { return velocity_; }

    /*!
     * \brief Return the maximum velocity vector \f$\mathrm{[m/s]}\f$ of the wall segment.
     */
    DimVector velocityMaximum() const
    { return velocityMaximum_; }

    /*!
     * \brief Return the velocity gradients \f$\mathrm{[1/s]}\f$ at the control volume center.
     */
    DimMatrix velocityGradients() const
    { return velocityGradients_; }

    /*!
     * \brief Return the wall distance \f$\mathrm{[m]}\f$ of the control volume.
     */
    Scalar wallDistance() const
    { return wallDistance_; }

    /*!
     * \brief Return the Karman constant
     */
    Scalar karmanConstant() const
    { return karmanConstant_; }

    /*!
     * \brief Return the wall friction velocity \f$\mathrm{[m/s]}\f$
     */
    Scalar uStar() const
    { return uStar_; }

    /*!
     * \brief Return the dimensionless wall distance \f$\mathrm{[-]}\f$.
     */
    Scalar yPlus() const
    { return yPlus_; }

    /*!
     * \brief Return the dimensionless velocity \f$\mathrm{[-]}\f$.
     */
    Scalar uPlus() const
    { return uPlus_; }

    /*!
     * \brief Return the dynamic eddy viscosity \f$\mathrm{[Pa s]}\f$ of the flow within the
     *        control volume.
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
     * \brief Return the kinematic eddy viscosity \f$\mathrm{[m^2/s]}\f$ of the flow within the
     *        control volume.
     */
    Scalar kinematicEddyViscosity() const
    { return dynamicEddyViscosity() / NavierStokesParentType::density(); }

    /*!
     * \brief Return the kinematic viscosity \f$\mathrm{[m^2/s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar kinematicViscosity() const
    { return NavierStokesParentType::viscosity() / NavierStokesParentType::density(); }

    /*!
     * \brief Calculates the eddy diffusivity \f$\mathrm{[m^2/s]}\f$ based
     *        on the kinematic eddy viscosity and the turbulent Schmidt number
     */
    template<class Problem>
    void calculateEddyDiffusivity(const Problem& problem)
    {
        eddyDiffusivity_ = kinematicEddyViscosity()
                           / problem.turbulentSchmidtNumber();
    }

    /*!
     * \brief Calculates the eddy thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ based
     *        on the kinematic eddy viscosity and the turbulent Prandtl number
     */
    template<class Problem, bool eB = enableEnergyBalance, typename std::enable_if_t<eB, int> = 0>
    void calculateEddyThermalConductivity(const Problem& problem)
    {
        eddyThermalConductivity_ = kinematicEddyViscosity()
                                   * NavierStokesParentType::density()
                                   * NavierStokesParentType::heatCapacity()
                                   / problem.turbulentPrandtlNumber();
    }

    //! \brief Eddy thermal conductivity is zero for isothermal model
    template<class Problem, bool eB = enableEnergyBalance, typename std::enable_if_t<!eB, int> = 0>
    void calculateEddyThermalConductivity(const Problem& problem)
    { eddyThermalConductivity_ = 0.0; }

    /*!
     * \brief Returns the eddy diffusivity \f$\mathrm{[m^2/s]}\f$
     */
    Scalar eddyDiffusivity() const
    { return eddyDiffusivity_; }

    /*!
     * \brief Returns the eddy thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
     */
    Scalar eddyThermalConductivity() const
    { return eddyThermalConductivity_; }

    /*!
    * \brief Returns the effective diffusion coefficient \f$\mathrm{[m^2/s]}\f$
    *
    * \param compIIdx the index of the component which diffusive
    * \param compJIdx the index of the component with respect to which compIIdx diffuses
    */
    Scalar effectiveDiffusivity(int compIIdx, int compJIdx) const
    {
        return NavierStokesParentType::diffusionCoefficient(compIIdx, compJIdx) + eddyDiffusivity();
    }

    /*!
     * \brief Returns the effective thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
     *        of the fluid-flow in the sub-control volume.
     */
    template<bool eB = enableEnergyBalance, typename std::enable_if_t<eB, int> = 0>
    Scalar effectiveThermalConductivity() const
    {
        return NavierStokesParentType::thermalConductivity() + eddyThermalConductivity();
    }

protected:
    /*!
     * \brief Sets the dynamic eddy viscosity \f$\mathrm{[Pa s]}\f$
     */
    Scalar setDynamicEddyViscosity_(Scalar value)
    { return dynamicEddyViscosity_  = value; }

    DimVector velocity_;
    DimVector velocityMaximum_;
    DimMatrix velocityGradients_;
    std::size_t elementID_;
    std::size_t wallElementID_;
    Scalar wallDistance_;
    Scalar karmanConstant_;
    Scalar uStar_;
    Scalar yPlus_;
    Scalar uPlus_;
    Scalar dynamicEddyViscosity_ = 0.0;
    Scalar eddyDiffusivity_ = 0.0;
    Scalar eddyThermalConductivity_ = 0.0;
};
}

#endif
