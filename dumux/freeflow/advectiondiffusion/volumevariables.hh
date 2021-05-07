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
 * \ingroup AdvectionDiffusionModel
 */
#ifndef DUMUX_ADVECTION_DIFFUSION_VOLUMEVARIABLES_HH
#define DUMUX_ADVECTION_DIFFUSION_VOLUMEVARIABLES_HH

#include <array>
#include <dune/common/exceptions.hh>
#include <dumux/freeflow/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowNCModel
 * \brief Volume variables for the single-phase, multi-component free-flow model.
 */
template <class Traits>
class AdvectionDiffusionVolumeVariables : public FreeFlowVolumeVariables< Traits, AdvectionDiffusionVolumeVariables<Traits> >
{
    using ThisType = AdvectionDiffusionVolumeVariables<Traits>;
    using ParentType = FreeFlowVolumeVariables<Traits, ThisType>;

    using Scalar = typename Traits::PrimaryVariables::value_type;

    static constexpr bool useMoles = Traits::ModelTraits::useMoles();
    static constexpr int numFluidComps = ParentType::numFluidComponents();
    using DiffusionCoefficients = typename Traits::DiffusionType::DiffusionCoefficientsContainer;

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
        // update parent type sets primary variables
        ParentType::update(elemSol, problem, element, scv);

        // the spatial params special to the tracer model
        fluidDensity_ = problem.spatialParams().fluidDensity(element, scv);
        fluidMolarMass_ = problem.spatialParams().fluidMolarMass(element, scv);

        for (int compIdx = 0; compIdx < ParentType::numFluidComponents(); ++compIdx)
        {
            moleOrMassFraction_[compIdx] = elemSol[0][Indices::componentIdx + compIdx];
            diffCoefficient_[compIdx] = FluidSystem::binaryDiffusionCoefficient(compIdx, problem, element, scv);
        }
    }

    /*!
     * \brief Returns the density \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx = 0) const
    { return fluidDensity_; }

    /*!
     * \brief Returns the average molar mass \f$\mathrm{[kg/mol]}\f$ of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar averageMolarMass(int phaseIdx = 0) const
    { return fluidMolarMass_; }

    /*!
     * \brief Returns the molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx = 0) const
    { return fluidDensity_ / fluidMolarMass_; }

    /*!
     * \brief Returns the mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    { return useMoles ? moleOrMassFraction_[compIdx]
                      : moleOrMassFraction_[compIdx] / FluidSystem::molarMass(compIdx) * fluidMolarMass_; }

    /*!
     * \brief Returns the mass fraction \f$\mathrm{[kg/kg]}\f$ of a component in the phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    { return useMoles ? moleOrMassFraction_[compIdx] * FluidSystem::molarMass(compIdx) / fluidMolarMass_
                      : moleOrMassFraction_[compIdx]; }

    /*!
     * \brief Returns the concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return moleFraction(phaseIdx, compIdx)*molarDensity(); }

    /*!
     * \brief Returns the binary diffusion coefficients for a phase in \f$[m^2/s]\f$.
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIIdx, int compJIdx) const
    { return diffCoefficient_[compJIdx]; }

    /*!
     * \brief Returns the effective diffusion coefficients for a phase in \f$[m^2/s]\f$.
     */
    Scalar effectiveDiffusionCoefficient(int phaseIdx, int compIIdx, int compJIdx) const
    { return diffusionCoefficient(0, compIIdx, compJIdx); }

protected:
    Scalar fluidDensity_, fluidMolarMass_;

    std::array<Scalar, numFluidComps> moleOrMassFraction_;
    // Binary diffusion coefficient
    std::array<Scalar, numFluidComps> diffCoefficient_;
};

} // end namespace Dumux

#endif
