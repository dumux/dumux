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
 * \ingroup TracerModel
 * \brief Quantities required by the tracer model in a control volume.
 */
#ifndef DUMUX_TRACER_VOLUME_VARIABLES_HH
#define DUMUX_TRACER_VOLUME_VARIABLES_HH

#include <cassert>
#include <array>
#include <type_traits>

#include <dune/common/std/type_traits.hh>

#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>

namespace Dumux {

namespace Detail {
// helper structs and functions detecting if the user-defined spatial params class
// has user-specified functions saturation() for multi-phase tracer.
template <typename T, typename ...Ts>
using SaturationDetector = decltype(std::declval<T>().spatialParams().saturation(std::declval<Ts>()...));

template<class T, typename ...Args>
static constexpr bool hasSaturation()
{ return Dune::Std::is_detected<SaturationDetector, T, Args...>::value; }

} // end namespace Detail

/*!
 * \ingroup TracerModel
 * \brief Contains the quantities which are constant within a
 *        finite volume for the tracer model.
 */
template <class Traits>
class TracerVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    static constexpr bool useMoles = Traits::ModelTraits::useMoles();
    using EffDiffModel = typename Traits::EffectiveDiffusivityModel;
    static constexpr int numFluidComps = ParentType::numFluidComponents();

public:
    //! Export the fluid system type
    using FluidSystem = typename Traits::FluidSystem;
    //! Export the solid state type
    using SolidState = typename Traits::SolidState;
    //! Export the indices
    using Indices = typename Traits::ModelTraits::Indices;

    /*!
     * \brief Updates all quantities for a given control volume.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv &scv)
    {
        // update parent type sets primary variables
        ParentType::update(elemSol, problem, element, scv);

        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);

        // the spatial params special to the tracer model
        fluidDensity_ = problem.spatialParams().fluidDensity(element, scv);
        fluidMolarMass_ = problem.spatialParams().fluidMolarMass(element, scv);

        if constexpr (Detail::hasSaturation<Problem, Element, Scv>())
            fluidSaturation_ = problem.spatialParams().saturation(element, scv);

        for (int compIdx = 0; compIdx < ParentType::numFluidComponents(); ++compIdx)
        {
            moleOrMassFraction_[compIdx] = this->priVars()[compIdx];
            diffCoeff_[compIdx] = FluidSystem::binaryDiffusionCoefficient(compIdx, problem, element, scv);
            effectiveDiffCoeff_[compIdx] = EffDiffModel::effectiveDiffusionCoefficient(*this, 0, 0, compIdx);
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
     * \brief Returns the phase state for the control volume.
     */
    const SolidState &solidState() const
    { return solidState_; }

    /*!
     * \brief Returns the saturation.
     *
     * This method is here for compatibility reasons with other models. The saturation
     * is always 1.0 in a one-phasic context, if two-phases or richards are considered,
     * the spatialParams serve as way to pass the saturation from the main-file to the
     * volVars and then to the localresidual for the tracer model.

     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx = 0) const
    { return fluidSaturation_ ; }

    /*!
     * \brief Returns the mobility.
     *
     * This method is here for compatibility reasons with other models. The mobility is always 1
     * for one-phasic models where the velocity field is given
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx = 0) const
    { return 1.0; }

    /*!
     * \brief Returns the molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx = 0) const
    { return fluidDensity_/fluidMolarMass_; }

    /*!
     * \brief Returns the mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    { return useMoles ? moleOrMassFraction_[compIdx] : moleOrMassFraction_[compIdx]/FluidSystem::molarMass(compIdx)*fluidMolarMass_; }

    /*!
     * \brief Returns the mass fraction \f$\mathrm{[kg/kg]}\f$ of a component in the phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    { return useMoles ? moleOrMassFraction_[compIdx]*FluidSystem::molarMass(compIdx)/fluidMolarMass_ : moleOrMassFraction_[compIdx]; }

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
    {
        if (phaseIdx != compIIdx) std::swap(compIIdx, compJIdx);
        assert(phaseIdx == 0);
        assert(phaseIdx == compIIdx);
        return diffCoeff_[compJIdx]; }

    /*!
     * \brief Returns the effective diffusion coefficients for a phase in \f$[m^2/s]\f$.
     */
    Scalar effectiveDiffusionCoefficient(int phaseIdx, int compIIdx, int compJIdx) const
    {
        if (phaseIdx != compIIdx) std::swap(compIIdx, compJIdx);
        assert(phaseIdx == 0);
        assert(phaseIdx == compIIdx);
        return effectiveDiffCoeff_[compJIdx];
    }

    /*!
     * \brief Return the average porosity \f$\mathrm{[-]}\f$ within the control volume.
     */
    Scalar porosity() const
    { return solidState_.porosity(); }

protected:
    SolidState solidState_;
    Scalar fluidDensity_, fluidMolarMass_;
    Scalar fluidSaturation_ = 1.0;

    std::array<Scalar, numFluidComps> diffCoeff_;
    std::array<Scalar, numFluidComps> effectiveDiffCoeff_;
    std::array<Scalar, numFluidComps> moleOrMassFraction_;
};

} // end namespace Dumux

#endif
