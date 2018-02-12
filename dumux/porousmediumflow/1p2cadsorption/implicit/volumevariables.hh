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
#ifndef DUMUX_1P2C_ADSORPTION_VOLUME_VARIABLES_HH
#define DUMUX_1P2C_ADSORPTION_VOLUME_VARIABLES_HH

#include <dumux/implicit/volumevariables.hh>

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup OnePTwoCAdsorptionModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are constant within a
 *        finite volume in the single-phase, two-component model.
 */
template <class TypeTag>
class OnePTwoCAdsorptionVolumeVariables : public ImplicitVolumeVariables<TypeTag>
{
    typedef ImplicitVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        phaseIdx        = Indices::phaseIdx,
        nPhaseIdx       = phaseIdx,
        phaseCompIdx    = Indices::phaseCompIdx,
        transportCompIdx= Indices::transportCompIdx,
        numComponents   = FluidSystem::numComponents,
        CH4Idx          = FluidSystem::CH4Idx,
        TCIdx           = FluidSystem::CO2Idx,

    };
    //indices of primary variables
    enum{
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar,dim> DimVector;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:

    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx,
                const bool isOldSol)
    {
        ParentType::update(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
        //calculate all secondary variables from the primary variables and store results in fluidstate
        completeFluidState(priVars, problem, element, fvGeometry, scvIdx, fluidState_);

        porosity_ = problem.spatialParams().porosity(element, fvGeometry, scvIdx);

        dispersivity_ = problem.spatialParams().dispersivity(element, fvGeometry, scvIdx);

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

        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(priVars, problem, element, fvGeometry, scvIdx, isOldSol);

        //Adsorption relevant things
        if (!GET_PARAM_FROM_GROUP(TypeTag, bool, SorptionCoefficients, useAdsorption))
                          {   // not using Adsorption (Standard if not otherwise mentioned in input
                              // set "useAdsorption = 0" in input file under [SorptionCoefficients] to enable Adsorption or delete from input file
                              adsorption_[numComponents] = {};
                          }
        else
        {   // using Adsorption
            //set "useAdsorption = 1" in input file under [SorptionCoefficients] to enable Adsorption
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                if (compIdx == CH4Idx)
                {
                    //check if Extended Langmuir formulation should be used
                    if (GET_PARAM_FROM_GROUP(TypeTag, bool, SorptionCoefficients, useEL))
                    {
                        adsorption_[compIdx] = adsorptionExtendedLangmuir(compIdx, problem);
                    }

                    //check if Simple Langmuir formulation should be used
                    else if (GET_PARAM_FROM_GROUP(TypeTag, bool, SorptionCoefficients, useL))
                    {
                        adsorption_[compIdx] = adsorptionSimpleLangmuir(compIdx, problem);
                    }

                    //check if Freundlich formulation should be used
                    else if (GET_PARAM_FROM_GROUP(TypeTag, bool, SorptionCoefficients, useF))
                    {
                        adsorption_[compIdx] = adsorptionFreundlich(compIdx, problem);
                    }

                    //check if BET formulation should be used
                    else if (GET_PARAM_FROM_GROUP(TypeTag, bool, SorptionCoefficients, useBET))
                    {
                        adsorption_[compIdx] = adsorptionBET(compIdx, problem);
                    }
                }
                else if (compIdx == TCIdx)
                {
                    //check if Extended Langmuir formulation should be used
                    if (GET_PARAM_FROM_GROUP(TypeTag, bool, SorptionCoefficients, useEL))
                    {
                        adsorption_[compIdx] = adsorptionExtendedLangmuir(compIdx, problem);
                    }

                    //check if Simple Langmuir formulation should be used
                    else if (GET_PARAM_FROM_GROUP(TypeTag, bool, SorptionCoefficients, useL))
                    {
                        adsorption_[compIdx] = adsorptionSimpleLangmuir(compIdx, problem);
                    }

                    //check if Freundlich formulation should be used
                    else if (GET_PARAM_FROM_GROUP(TypeTag, bool, SorptionCoefficients, useF))
                    {
                        adsorption_[compIdx] = adsorptionFreundlich(compIdx, problem);
                    }

                    //check if BET formulation should be used
                    else if (GET_PARAM_FROM_GROUP(TypeTag, bool, SorptionCoefficients, useBET))
                    {
                        adsorption_[compIdx] = adsorptionBET(compIdx, problem);
                    }
                }
                else
                    adsorption_[compIdx] = 0;
            }
        }
    }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     */
    static void completeFluidState(const PrimaryVariables& priVars,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const int scvIdx,
                                   FluidState& fluidState)
    {
        Scalar t = Implementation::temperature_(priVars, problem, element,
                                                fvGeometry, scvIdx);
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
        Scalar h = Implementation::enthalpy_(fluidState, paramCache, phaseIdx);
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
     */
    Scalar density() const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     */
    Scalar molarDensity() const
    { return fluidState_.molarDensity(phaseIdx);}

    /*!
     * \brief Return mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
     *
     * \param compIdx The index of the component
     */
    Scalar moleFraction(int compIdx) const
    { return fluidState_.moleFraction(phaseIdx, (compIdx==0)?phaseCompIdx:transportCompIdx); }

    /*!
     * \brief Return mass fraction \f$\mathrm{[kg/kg]}\f$ of a component in the phase.
     *
     * \param compIdx The index of the component
     */
    Scalar massFraction(int compIdx) const
    { return fluidState_.massFraction(phaseIdx, (compIdx==0)?phaseCompIdx:transportCompIdx); }

    /*!
     * \brief Return concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param compIdx The index of the component
     */
    Scalar molarity(int compIdx) const
    { return fluidState_.molarity(phaseIdx, (compIdx==0)?phaseCompIdx:transportCompIdx); }

    /*!
     * \brief Return the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     */
    Scalar pressure() const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ in the fluid.
     */
    Scalar diffCoeff() const
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
    { return fluidState_.temperature(phaseIdx); }

    /*!
     * \brief Return the dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of a given phase
     *        within the control volume.
     */
    Scalar viscosity() const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Return the average porosity \f$\mathrm{[-]}\f$ within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the adsorbed volume \f$\mathrm{a_ads}\f$ in \f$\mathrm{[mol/m^3]}\f$ of a component,
     * which is then added to the n-Phase storage term.
     *
     * \param compIdx The component index
     */
    Scalar adsorption(int compIdx) const
    {
        return adsorption_[compIdx];
    }

    /*!
     * \brief Returns the amount of adsorbate in the adsorbent at equilibrium.
     *      using Extended Langmuir Isotherm
     *
     * \param compIdx The component index
     */
    Scalar adsorptionExtendedLangmuir(int compIdx, const Problem &problem)
    {
        // using Extended Langmuir (not Gibbs conform!), which is the standard if not ohterwise mentioned in input
        // "useEL=1", "useL=0", "useF=0","useBET=0" in input file under [SorptionCoefficients] to enable Extended Langmuir or delete line
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            if (compIdx == CH4Idx)
            {
                adsorptionExtendedLangmuir_[compIdx] = problem.V(compIdx) * problem.b(compIdx) *
                    this->fluidState_.moleFraction(nPhaseIdx, compIdx)* this->fluidState_.pressure(nPhaseIdx)  /
                    (1 + problem.b(compIdx) * this->fluidState_.moleFraction(nPhaseIdx, compIdx)*
                    this->fluidState_.pressure(nPhaseIdx) + problem.b(TCIdx) * this->fluidState_.moleFraction(nPhaseIdx, TCIdx) *
                    this->fluidState_.pressure(nPhaseIdx));
            }
            else if (compIdx == TCIdx)// Extended Langmuir (not Gibbs conform!)
            {
                adsorptionExtendedLangmuir_[compIdx] = problem.V(compIdx) * problem.b(compIdx) *
                    this->fluidState_.moleFraction(nPhaseIdx, compIdx) * this->fluidState_.pressure(nPhaseIdx) /
                    (1 + problem.b(compIdx) * this->fluidState_.moleFraction(nPhaseIdx, compIdx) *
                    this->fluidState_.pressure(nPhaseIdx) + problem.b(CH4Idx)* this->fluidState_.moleFraction(nPhaseIdx, CH4Idx) *
                    this->fluidState_.pressure(nPhaseIdx));
            }
            else
            {
                adsorptionExtendedLangmuir_[compIdx] = 0;
            }
        }
        return adsorptionExtendedLangmuir_[compIdx];
    }

    /*!
     * \brief Returns the amount of adsorbate in the adsorbent at equilibrium.
     *      using Single Langmuir Isotherm
     *
     * \param compIdx The component index
     */
    Scalar adsorptionSimpleLangmuir(int compIdx, const Problem &problem)
    {
        // using simple Langmuir Adsorption and not using EL formulation now
        // set "useEL=0", "useL=1", "useF=0", "useBET=0" in input file under [SorptionCoefficients] to use Simple Langmuir Adsorption

        //Standard Langmuir Adsorption formulation
        adsorptionSimpleLangmuir_[compIdx] = problem.V(compIdx) * problem.b(compIdx) *
            this->fluidState_.moleFraction(nPhaseIdx, compIdx) * this->fluidState_.pressure(nPhaseIdx) /
            (1 + problem.b(compIdx) * this->fluidState_.moleFraction(nPhaseIdx, compIdx) *
            this->fluidState_.pressure(nPhaseIdx) );

        return adsorptionSimpleLangmuir_[compIdx];
    }

    /*!
     * \brief Returns the amount of adsorbate in the adsorbent at equilibrium.
     *      using Freundlich Isotherm
     *
     * \param compIdx The component index
     */
    Scalar adsorptionFreundlich(int compIdx, const Problem &problem)
    {
        //using Freundlich Adsorption
        //set "useEL=0", "useL=0", "useF=1", "useBET=0" in input file under [SorptionCoefficients] to use Freundlich Adsorption

        //Freundlich Adsorption Formulation
        adsorptionFreundlich_[compIdx] = problem.kF() * std::pow((this->fluidState_.moleFraction(nPhaseIdx, compIdx) *
            this->fluidState_.pressure(nPhaseIdx)), (1/problem.nF()));

        return adsorptionFreundlich_[compIdx];
    }

     /*!
     * \brief Returns the amount of adsorbate in the adsorbent at equilibrium.
     *      using BET Isotherm
     *
     * \param compIdx The component index
     */
    Scalar adsorptionBET(int compIdx, const Problem &problem)
    {
        //using BET Adsorption
        //set "useEL=0", "useL=0", "useF=0", "useBET=1" in input file under [Sorption Coefficients] to use Freundlich Adsorption

        //BET Adsorption Formulation
        adsorptionBET_[compIdx] = problem.qsatBET() * problem.cBET() * this->fluidState_.moleFraction(nPhaseIdx, compIdx) *
            this->fluidState_.pressure(nPhaseIdx) / ((problem.csatBET() - this->fluidState_.moleFraction(nPhaseIdx, compIdx) *
            this->fluidState_.pressure(nPhaseIdx)) * (1 + (problem.cBET() - 1) * (this->fluidState_.moleFraction(nPhaseIdx, compIdx) *
            this->fluidState_.pressure(nPhaseIdx) / problem.csatBET() )));

        return adsorptionBET_[compIdx];
    }

protected:
    static Scalar temperature_(const PrimaryVariables &priVars,
                               const Problem& problem,
                               const Element &element,
                               const FVElementGeometry &fvGeometry,
                               const int scvIdx)
    {
        return problem.temperatureAtPos(fvGeometry.subContVol[scvIdx].global);
    }

    template<class ParameterCache>
    static Scalar enthalpy_(const FluidState& fluidState,
                            const ParameterCache& paramCache,
                            const int phaseIdx)
    {
        return 0;
    }

    /*!
     * \brief Called by update() to compute the energy related quantities.
     */
    void updateEnergy_(const PrimaryVariables &priVars,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &fvGeometry,
                       const int scvIdx,
                       const bool isOldSol)
    { }

    Scalar porosity_;    //!< Effective porosity within the control volume
    GlobalPosition dispersivity_;
    Scalar diffCoeff_;
    FluidState fluidState_;
    Scalar adsorption_[numComponents] = {};
    Scalar adsorptionExtendedLangmuir_[numComponents] = {};
    Scalar adsorptionSimpleLangmuir_[numComponents] = {};
    Scalar adsorptionFreundlich_[numComponents] = {};
    Scalar adsorptionBET_[numComponents] = {};
private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}// end namespace

#endif
