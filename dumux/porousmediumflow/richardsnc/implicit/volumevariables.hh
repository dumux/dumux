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
 * \brief @copybrief Dumux::RichardsNCVolumeVariables
 */
#ifndef DUMUX_RICHARDSNC_VOLUME_VARIABLES_HH
#define DUMUX_RICHARDSNC_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/richards/implicit/volumevariables.hh>

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup RichardsNCModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are constant within a
 *        finite volume in the Richards, n-component model.
 */
template <class TypeTag>
class RichardsNCVolumeVariables : public RichardsVolumeVariables<TypeTag>
{
    using ParentType = RichardsVolumeVariables<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using PermeabilityType = typename SpatialParams::PermeabilityType;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pressureIdx = Indices::pressureIdx
    };

    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    static const int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        //calculate all secondary variables from the primary variables and store results in fluidstate
        Implementation::completeFluidState(elemSol, problem, element, scv, this->fluidState_);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(this->fluidState_, wPhaseIdx);

        const int compIIdx = wPhaseIdx;
        for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            if(compIIdx != compJIdx)
                setDiffusionCoefficient_(compJIdx,
                                         FluidSystem::binaryDiffusionCoefficient(this->fluidState_,
                                                                                 paramCache,
                                                                                 wPhaseIdx,
                                                                                 compIIdx,
                                                                                 compJIdx));
    }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     */
    static void completeFluidState(const ElementSolutionVector &elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume &scv,
                                   FluidState& fluidState)
    {
        ParentType::completeFluidState(elemSol, problem, element, scv, fluidState);

        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv);

        // set the mole/mass fractions
        if(useMoles)
        {
            Scalar sumSecondaryFractions = 0.0;
            for (int compIdx = 1; compIdx < numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(wPhaseIdx, compIdx, priVars[compIdx]);
                sumSecondaryFractions += priVars[compIdx];
            }
            fluidState.setMoleFraction(wPhaseIdx, 0, 1.0 - sumSecondaryFractions);
        }
        else
        {
            for (int compIdx = 1; compIdx < numComponents; ++compIdx)
                fluidState.setMassFraction(wPhaseIdx, compIdx, priVars[compIdx]);
        }
    }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar molarDensity(const int phaseIdx = wPhaseIdx) const
    { return phaseIdx == wPhaseIdx ? this->fluidState_.molarDensity(phaseIdx) : 0.0; }

    /*!
     * \brief Return mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
     *
     * \param compIdx The index of the component
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar moleFraction(const int phaseIdx, const int compIdx) const
    { return phaseIdx == wPhaseIdx ? this->fluidState_.moleFraction(phaseIdx, compIdx) : 0.0; }

    /*!
     * \brief Return mass fraction \f$\mathrm{[kg/kg]}\f$ of a component in the phase.
     *
     * \param compIdx The index of the component
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar massFraction(const int phaseIdx, const int compIdx) const
    { return phaseIdx == wPhaseIdx ? this->fluidState_.massFraction(phaseIdx, compIdx) : 0.0; }

    /*!
     * \brief Return concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param compIdx The index of the component
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar molarity(const int phaseIdx, const int compIdx) const
    { return phaseIdx == wPhaseIdx ? this->fluidState_.molarity(phaseIdx, compIdx) : 0.0; }

    /*!
     * \brief Return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ in the fluid.
     */
    Scalar diffusionCoefficient(const int phaseIdx, const int compIdx) const
    {
        assert(phaseIdx == wPhaseIdx);
        assert(compIdx > wPhaseIdx);
        return diffCoefficient_[compIdx-1];
    }

    /*!
     * \brief Returns the dispersivity of the fluid's streamlines.
     */
    const GlobalPosition &dispersivity() const
    { return dispersivity_; }

private:
    void setDiffusionCoefficient_(int compIdx, Scalar d)
    {
        assert(compIdx > wPhaseIdx);
        diffCoefficient_[compIdx-1] = d;
    }

    std::array<Scalar, numComponents-1> diffCoefficient_;
    GlobalPosition dispersivity_;
};

} // end namespace Dumux

#endif
