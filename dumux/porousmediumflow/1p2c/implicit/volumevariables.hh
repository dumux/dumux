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
 * \brief Quantities required by the single-phase, two-component box
 *        model defined on a vertex.
 */
#ifndef DUMUX_1P2C_VOLUME_VARIABLES_HH
#define DUMUX_1P2C_VOLUME_VARIABLES_HH

#include <dumux/discretization/volumevariables.hh>

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup OnePTwoCModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are constant within a
 *        finite volume in the single-phase, two-component model.
 *
 * \note The return functions for the fluid state variables always forward to the actual
 *       fluid state using the phaseIdx from the DuMuX property system. Furthermore, the
 *       default value is not used, but is only here to enable calling these functions
 *       without handing in a phase index (as in a single-phasic context there is only one phase).
 *       This way one can use two-phase fluid systems for this one-phasic flow and transport
 *       model by specifying which phase is present through the DuMuX property system.
 */
template <class TypeTag>
class OnePTwoCVolumeVariables : public ImplicitVolumeVariables<TypeTag>
{
    using ParentType = ImplicitVolumeVariables<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum
    {
        phaseIdx = Indices::phaseIdx,
        phaseCompIdx = Indices::phaseCompIdx,
        transportCompIdx = Indices::transportCompIdx,

        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx
    };

    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimVector = Dune::FieldVector<Scalar,dim>;
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

public:

    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    {
        ParentType::update(priVars, problem, element, scv);

        //calculate all secondary variables from the primary variables and store results in fluidstate
        completeFluidState(priVars, problem, element, scv, fluidState_);

        porosity_ = problem.spatialParams().porosity(scv);
        dispersivity_ = problem.spatialParams().dispersivity(element, scv);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState_, phaseIdx);

        diffCoeff_ = FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                             paramCache,
                                                             phaseIdx,
                                                             phaseCompIdx,
                                                             transportCompIdx);

        Valgrind::CheckDefined(porosity_);
        Valgrind::CheckDefined(dispersivity_);
        Valgrind::CheckDefined(diffCoeff_);
    }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     */
    static void completeFluidState(const PrimaryVariables& priVars,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume &scv,
                                   FluidState& fluidState)
    {
        Scalar t = ParentType::temperature(priVars, problem, element, scv);
        fluidState.setTemperature(t);
        fluidState.setSaturation(phaseIdx, 1.);

        fluidState.setPressure(phaseIdx, priVars[pressureIdx]);

        if(useMoles)
        {
            fluidState.setMoleFraction(phaseIdx, phaseCompIdx, 1 - priVars[massOrMoleFracIdx]);
            fluidState.setMoleFraction(phaseIdx, transportCompIdx, priVars[massOrMoleFracIdx]);
        }
        else
        {
            // setMassFraction() has only to be called 1-numComponents times
            fluidState.setMassFraction(phaseIdx, transportCompIdx, priVars[massOrMoleFracIdx]);
        }

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, phaseIdx);

        Scalar value;
        value = FluidSystem::density(fluidState, paramCache, phaseIdx);
        fluidState.setDensity(phaseIdx, value);
        value = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
        fluidState.setViscosity(phaseIdx, value);

        // compute and set the enthalpy
        Scalar h = Implementation::enthalpy(fluidState, paramCache, phaseIdx);
        fluidState.setEnthalpy(phaseIdx, h);
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar density(int pIdx = 0) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Return the saturation
     *
     * This method is here for compatibility reasons with other models. The saturation
     * is always 1.0 in a one-phasic context.
     */
    Scalar saturation(int pIdx = 0) const
    { return 1.0; }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar molarDensity(int pIdx = 0) const
    { return fluidState_.molarDensity(phaseIdx); }

    /*!
     * \brief Return mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
     *
     * \param compIdx The index of the component
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar moleFraction(int pIdx, int compIdx) const
    { return fluidState_.moleFraction(phaseIdx, (compIdx==0)?phaseCompIdx:transportCompIdx); }

    /*!
     * \brief Return mass fraction \f$\mathrm{[kg/kg]}\f$ of a component in the phase.
     *
     * \param compIdx The index of the component
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar massFraction(int pIdx, int compIdx) const
    { return fluidState_.massFraction(phaseIdx, (compIdx==0)?phaseCompIdx:transportCompIdx); }

    /*!
     * \brief Return concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param compIdx The index of the component
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar molarity(int pIdx, int compIdx) const
    { return fluidState_.molarity(phaseIdx, (compIdx==0)?phaseCompIdx:transportCompIdx); }

    /*!
     * \brief Return the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar pressure(int pIdx = 0) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ in the fluid.
     */
    Scalar diffusionCoefficient(int pIdx, int compIdx) const
    { return diffCoeff_; }

    /*!
     * \brief Returns the dispersivity of the fluid's streamlines.
     */
    const GlobalPosition &dispersivity() const
    { return dispersivity_; }

    /*!
     * \brief Return temperature \f$\mathrm{[K]}\f$ inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

    /*!
     * \brief Returns the mobility \f$\mathrm{[1/(Pa s)]}\f$.
     *
     * The term mobility is usually not employed in the one phase context.
     * The method is here for compatibility reasons with other models.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar mobility(int pIdx = 0) const
    { return 1.0/fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Return the dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of a given phase
     *        within the control volume.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar viscosity(int pIdx = 0) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Return the average porosity \f$\mathrm{[-]}\f$ within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

protected:
    Scalar porosity_;    //!< Effective porosity within the control volume
    GlobalPosition dispersivity_;
    Scalar diffCoeff_;
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}// end namespace

#endif
