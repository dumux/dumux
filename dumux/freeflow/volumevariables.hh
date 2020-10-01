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
 * \ingroup FreeflowModels
 *
 * \copydoc Dumux::FreeFlowVolumeVariables
 */
#ifndef DUMUX_FREEFLOW_VOLUME_VARIABLES_HH
#define DUMUX_FREEFLOW_VOLUME_VARIABLES_HH

#include <dumux/material/fluidstates/immiscible.hh>

namespace Dumux {

// forward declaration
template <class Traits, class Impl, bool enableEnergyBalance>
class FreeFlowVolumeVariablesImplementation;

/*!
 * \ingroup FreeflowModels
 * \brief Volume variables for free-flow models.
 *        The class is specialized for isothermal and non-isothermal models.
 */
template <class Traits, class Impl>
using FreeFlowVolumeVariables = FreeFlowVolumeVariablesImplementation<Traits, Impl, Traits::ModelTraits::enableEnergyBalance()>;

/*!
 * \ingroup FreeflowModels
 * \brief Volume variables for isothermal free-flow models.
 */
template <class Traits, class Impl>
class FreeFlowVolumeVariablesImplementation<Traits, Impl, false>
{
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;
    //! export the type encapsulating primary variable indices
    using Indices = typename Traits::ModelTraits::Indices;

    //! return number of phases considered by the model
    static constexpr int numFluidPhases() { return Traits::ModelTraits::numFluidPhases(); }
    //! return number of components considered by the model
    static constexpr int numFluidComponents() { return Traits::ModelTraits::numFluidComponents(); }

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
        priVars_ = elemSol[scv.localDofIndex()];
        extrusionFactor_ = problem.extrusionFactor(element, scv, elemSol);
    }

    /*!
     * \brief Return how much the sub-control volume is extruded.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    /*!
     * \brief Return a component of primary variable vector
     *
     * \param pvIdx The index of the primary variable of interest
     */
    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    //! The temperature is obtained from the problem as a constant for isothermal models
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    static Scalar temperature(const ElementSolution &elemSol,
                              const Problem& problem,
                              const Element &element,
                              const SubControlVolume &scv)
    {
        return problem.temperatureAtPos(scv.dofPosition());
    }

    //! The phase enthalpy is zero for isothermal models
    //! This is needed for completing the fluid state
    template<class FluidState, class ParameterCache>
    static Scalar enthalpy(const FluidState& fluidState,
                           const ParameterCache& paramCache)
    {
        return 0;
    }

protected:
    const Impl &asImp_() const { return *static_cast<const Impl*>(this); }
    Impl &asImp_() { return *static_cast<Impl*>(this); }
    PrimaryVariables priVars_;
    Scalar extrusionFactor_;
};

/*!
 * \ingroup FreeflowModels
 * \brief Volume variables for the non-isothermal free-flow models.
 */
template <class Traits, class Impl>
class FreeFlowVolumeVariablesImplementation<Traits, Impl, true>
: public FreeFlowVolumeVariablesImplementation<Traits, Impl, false>
{
    using ParentType = FreeFlowVolumeVariablesImplementation<Traits, Impl, false>;
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;
    //! export the underlying fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! export the type encapsulating primary variable indices
    using Indices = typename Traits::ModelTraits::Indices;


    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume which is inside the element
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv &scv)
    {
        ParentType::update(elemSol, problem, element, scv);
    }

    /*!
     * \brief Returns the total internal energy of a phase in the
     *        sub-control volume.
     */
    Scalar internalEnergy(int phaseIdx = 0) const
    { return ParentType::asImp_().fluidState().internalEnergy(0); }

    /*!
     * \brief Returns the total enthalpy of a phase in the sub-control
     *        volume.
     */
    Scalar enthalpy(int phaseIdx = 0) const
    { return ParentType::asImp_().fluidState().enthalpy(0); }

    /*!
     * \brief Returns the component enthalpy \f$\mathrm{[J/(kg*K)]}\f$ in the sub-control volume.
     */
    Scalar componentEnthalpy(unsigned int compIdx) const
    { return FluidSystem::componentEnthalpy(ParentType::asImp_().fluidState(), 0, compIdx); }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
     *        of the fluid phase in the sub-control volume.
     */
    Scalar thermalConductivity() const
    { return FluidSystem::thermalConductivity(ParentType::asImp_().fluidState(), 0); }

    /*!
     * \brief Returns the effective thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
     *        of the fluid-flow in the sub-control volume.
     */
    Scalar effectiveThermalConductivity() const
    {
        return thermalConductivity();
    }

    /*!
     * \brief Return the specific isobaric heat capacity \f$\mathrm{[J/(kg*K)]}\f$
     *        in the sub-control volume.
     */
    Scalar heatCapacity() const
    { return FluidSystem::heatCapacity(ParentType::asImp_().fluidState(), 0); }

    //! The temperature is a primary variable for non-isothermal models
    template<class ElemSol, class Problem, class Element, class Scv>
    static Scalar temperature(const ElemSol &elemSol,
                              const Problem& problem,
                              const Element &element,
                              const Scv &scv)
    {
        return elemSol[scv.localDofIndex()][Indices::temperatureIdx];
    }

    //! The phase enthalpy is zero for isothermal models
    //! This is needed for completing the fluid state
    template<class FluidState, class ParameterCache>
    static Scalar enthalpy(const FluidState& fluidState,
                           const ParameterCache& paramCache)
    {
        return FluidSystem::enthalpy(fluidState, paramCache, 0);
    }
};
} // end namespace Dumux

#endif
