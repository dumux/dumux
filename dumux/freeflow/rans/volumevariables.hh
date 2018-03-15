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

#include <dumux/common/properties.hh>
#include <dumux/material/fluidstates/immiscible.hh>

namespace Dumux
{

// forward declaration
template <class TypeTag, bool enableEnergyBalance>
class RANSVolumeVariablesImplementation;

/*!
 * \ingroup RANSModel
 * \brief Volume variables for the single-phase Reynolds-Averaged Navier-Stokes models.
 *        The class is specialized for isothermal and non-isothermal models.
 */
template <class TypeTag>
using RANSVolumeVariables = RANSVolumeVariablesImplementation<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

/*!
 * \ingroup RANSModel
 * \brief Volume variables for the isothermal single-phase Reynolds-Averaged Navier-Stokes models.
 */
template <class TypeTag>
class RANSVolumeVariablesImplementation<TypeTag, false>
: public NavierStokesVolumeVariablesImplementation<TypeTag, false>
{
    using ParentType = NavierStokesVolumeVariablesImplementation<TypeTag, false>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    enum { dimWorld = GridView::dimensionworld };
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    static const int defaultPhaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);

public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \brief Update all quantities for a given control volume
     *
     * Wall related quantities are stored and the calculateEddyViscosity(...)
     * function of the turbulence model implementation is called.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        using std::abs;
        using std::max;
        using std::sqrt;
        ParentType::update(elemSol, problem, element, scv);

        // calculate characteristic properties of the turbulent flow
        elementID_ = problem.fvGridGeometry().elementMapper().index(element);
        wallElementID_ = problem.wallElementIDs_[elementID_];
        wallDistance_ = problem.wallDistances_[elementID_];
        velocity_ = problem.velocity_[elementID_];
        velocityGradients_ = problem.velocityGradients_[elementID_];
        unsigned int flowNormalAxis = problem.flowNormalAxis_[elementID_];
        unsigned int wallNormalAxis = problem.wallNormalAxis_[elementID_];
        Scalar uStar = sqrt(problem.kinematicViscosity_[wallElementID_]
                            * abs(problem.velocityGradients_[wallElementID_][flowNormalAxis][wallNormalAxis]));
        yPlus_ = wallDistance_ * uStar / asImp_().kinematicViscosity();
        yPlus_ = max(yPlus_, 1e-10); // zero values lead to numerical problems in some turbulence models
        uPlus_ = velocity_[flowNormalAxis] / max(uStar, 1e-10);
    };

    /*!
     * \brief Return the velocity vector \f$\mathrm{[m/s]}\f$ at the control volume center.
     */
    DimVector velocity() const
    { return velocity_; }

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
     * \brief Set the values of the dynamic eddy viscosity \f$\mathrm{[Pa s]}\f$ within the
     *        control volume.
     */
    void setDynamicEddyViscosity(Scalar value)
    { dynamicEddyViscosity_ = value; }

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
    { return asImp_().viscosity() + dynamicEddyViscosity(); }

    /*!
     * \brief Return the kinematic eddy viscosity \f$\mathrm{[m^2/s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar kinematicViscosity() const
    { return asImp_().viscosity() / asImp_().density(); }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

protected:
    DimVector velocity_;
    DimMatrix velocityGradients_;
    Scalar dynamicEddyViscosity_;
    unsigned int elementID_;
    unsigned int wallElementID_;
    Scalar wallDistance_;
    Scalar yPlus_;
    Scalar uPlus_;
};

/*!
 * \ingroup RANSModel
 * \brief Volume variables for the non-isothermal single-phase Reynolds-Averaged Navier-Stokes models.
 */
template <class TypeTag>
class RANSVolumeVariablesImplementation<TypeTag, true>
: public NavierStokesVolumeVariablesImplementation<TypeTag, true>,
  public RANSVolumeVariablesImplementation<TypeTag, false>
{
    using ParentTypeNonIsothermal = NavierStokesVolumeVariablesImplementation<TypeTag, true>;
    using ParentTypeIsothermal = RANSVolumeVariablesImplementation<TypeTag, false>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    static const int defaultPhaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
//     static const int temperatureIdx = Indices::temperatureIdx;

public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    {
        ParentTypeIsothermal::update(elemSol, problem, element, scv);
        ParentTypeNonIsothermal::update(elemSol, problem, element, scv);
        // TODO: convert eddyViscosity to eddyThermalConductivity
    }

    /*!
     * \brief Returns the eddy thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
     *        of the flow phase in the sub-control volume.
     */
    Scalar eddyThermalConductivity() const
    { return eddyThermalConductivity_; }

    /*!
     * \brief Returns the effective thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
     *        of the fluid-flow in the sub-control volume.
     */
    Scalar effectiveThermalConductivity() const
    {
        return FluidSystem::thermalConductivity(this->fluidState_)
               + eddyThermalConductivity();
    }

protected:
    Scalar eddyThermalConductivity_;
};
}

#endif
