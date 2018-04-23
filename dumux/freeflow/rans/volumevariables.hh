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

// forward declaration
template <class Traits, class Impl, bool enableEnergyBalance>
class RANSVolumeVariablesImplementation;

/*!
 * \ingroup RANSModel
 * \brief Volume variables for the single-phase Reynolds-Averaged Navier-Stokes models.
 *        The class is specialized for isothermal and non-isothermal models.
 */
template <class Traits, class Impl>
using RANSVolumeVariables = RANSVolumeVariablesImplementation<Traits, Impl, Traits::ModelTraits::enableEnergyBalance()>;

/*!
 * \ingroup RANSModel
 * \brief Volume variables for the isothermal single-phase Reynolds-Averaged Navier-Stokes models.
 */
template <class Traits, class Impl>
class RANSVolumeVariablesImplementation<Traits, Impl, false>
{
    using Scalar = typename Traits::PrimaryVariables::value_type;

    enum { dimWorld = Traits::ModelTraits::dim() };
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:

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
     * \brief Return the kinematic eddy viscosity \f$\mathrm{[m^2/s]}\f$ of the flow within the
     *        control volume.
     */
    Scalar kinematicEddyViscosity() const
    { return asImp_().dynamicEddyViscosity() / asImp_().density(); }

    /*!
     * \brief Return the kinematic viscosity \f$\mathrm{[m^2/s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar kinematicViscosity() const
    { return asImp_().viscosity() / asImp_().density(); }

protected:
    DimVector velocity_;
    DimVector velocityMaximum_;
    DimMatrix velocityGradients_;
    std::size_t elementID_;
    std::size_t wallElementID_;
    Scalar wallDistance_;
    Scalar uStar_;
    Scalar yPlus_;
    Scalar uPlus_;

    const Impl &asImp_() const { return *static_cast<const Impl*>(this); }
    Impl &asImp_() { return *static_cast<Impl*>(this); }
};

/*!
 * \ingroup RANSModel
 * \brief Volume variables for the non-isothermal single-phase Reynolds-Averaged Navier-Stokes models.
 */
template <class Traits, class Impl>
class RANSVolumeVariablesImplementation<Traits, Impl, true>
: public RANSVolumeVariablesImplementation<Traits, Impl, false>
{
    using ParentType = RANSVolumeVariablesImplementation<Traits, Impl, false>;
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:

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
        ParentType::updateRANSProperties(elemSol, problem, element, scv);
        calculateEddyThermalConductivity(problem);
    }

    /*!
     * \brief Calculates the eddy thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ based
     *        on the kinematic eddy viscosity and the turbulent prandtl number
     */
    template<class Problem>
    void calculateEddyThermalConductivity(const Problem& problem)
    {
        static const auto turbulentPrandtlNumber
            = getParamFromGroup<Scalar>(problem.paramGroup(),
                                        "RANS.TurbulentPrandtlNumber", 1.0);

        eddyThermalConductivity_ = ParentType::kinematicEddyViscosity()
                                   * ParentType::asImp_().density()
                                   * ParentType::asImp_().heatCapacity()
                                   / turbulentPrandtlNumber;
    }

    /*!
     * \brief Returns the eddy thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
     *        of the flow in the sub-control volume.
     */
    Scalar eddyThermalConductivity() const
    { return eddyThermalConductivity_; }

protected:
    Scalar eddyThermalConductivity_;
};
}

#endif
