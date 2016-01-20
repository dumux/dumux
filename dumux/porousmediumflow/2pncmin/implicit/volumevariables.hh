// -**- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 *
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase, n-component mineralization model.
 */
#ifndef DUMUX_2PNCMin_VOLUME_VARIABLES_HH
#define DUMUX_2PNCMin_VOLUME_VARIABLES_HH

#include <dumux/implicit/model.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>
#include <dumux/common/math.hh>
#include <vector>
#include <iostream>

#include "properties.hh"
#include "indices.hh"
#include <dumux/material/constraintsolvers/computefromreferencephase2pncmin.hh>
#include <dumux/material/constraintsolvers/miscible2pnccomposition.hh>
#include <dumux/porousmediumflow/2pnc/implicit/volumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPNCMinModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, n-component model.
 */
template <class TypeTag>
class TwoPNCMinVolumeVariables : public TwoPNCVolumeVariables<TypeTag>
{
    typedef TwoPNCVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        dim = GridView::dimension,
        dimWorld=GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numSPhases =  GET_PROP_VALUE(TypeTag, NumSPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        numMajorComponents = GET_PROP_VALUE(TypeTag, NumMajorComponents),

        // formulations
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        plSg = TwoPNCFormulation::plSg,
        pgSl = TwoPNCFormulation::pgSl,

        // phase indices
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        // component indices
        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        // phase presence enums
        nPhaseOnly = Indices::nPhaseOnly,
        wPhaseOnly = Indices::wPhaseOnly,
        bothPhases = Indices::bothPhases,

        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        useSalinity = GET_PROP_VALUE(TypeTag, useSalinity)
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename Grid::ctype CoordScalar;
    typedef Dumux::Miscible2pNCComposition<Scalar, FluidSystem> Miscible2pNCComposition;
    typedef Dumux::ComputeFromReferencePhase2pNCMin<Scalar, FluidSystem> ComputeFromReferencePhase2pNCMin;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };
public:

    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx,
                bool isOldSol)
    {
        ParentType::update(priVars,
                           problem,
                           element,
                           fvGeometry,
                           scvIdx,
                           isOldSol);

        completeFluidState(priVars, problem, element, fvGeometry, scvIdx, this->fluidState_, isOldSol);

    /////////////
        // calculate the remaining quantities
        /////////////

    // porosity evaluation
    initialPorosity_ = problem.spatialParams().porosity(element, fvGeometry, scvIdx);
    minimumPorosity_ = problem.spatialParams().porosityMin(element, fvGeometry, scvIdx);


    sumPrecipitates_ = 0.0;
    for(int sPhaseIdx = 0; sPhaseIdx < numSPhases; ++sPhaseIdx)
    {
       precipitateVolumeFraction_[sPhaseIdx] = priVars[numComponents + sPhaseIdx];
       sumPrecipitates_+= precipitateVolumeFraction_[sPhaseIdx];
    }

//         for(int sPhaseIdx = 0; sPhaseIdx < numSPhases; ++sPhaseIdx)
//     {
//         Chemistry chemistry; // the non static functions can not be called without abject
//         saturationIdx_[sPhaseIdx] = chemistry.omega(sPhaseIdx);
//     }
// TODO/FIXME: The salt crust porosity is not clearly defined. However form literature review it is
//    found that the salt crust have porosity of approx. 10 %. Thus we restrict the decrease in porosity
//    to this limit. Moreover in the Problem files the precipitation should also be made dependent on local
//    porosity value, as the porous media media properties change related to salt precipitation will not be
//    accounted otherwise.

//      this->porosity_ = initialPorosity_ - sumPrecipitates_;

     this->porosity_ = std::max(minimumPorosity_, std::max(0.0, initialPorosity_ - sumPrecipitates_));

   salinity_= 0.0;
   moleFractionSalinity_ = 0.0;
   for (int compIdx = numMajorComponents; compIdx< numComponents; compIdx++)    //sum of the mass fraction of the components
   {
       if(this->fluidState_.moleFraction(wPhaseIdx, compIdx)> 0)
       {
          salinity_+= this->fluidState_.massFraction(wPhaseIdx, compIdx);
          moleFractionSalinity_ += this->fluidState_.moleFraction(wPhaseIdx, compIdx);
       }
    }

// TODO/FIXME: Different relations for the porosoty-permeability changes are given here. We have to fins a way
//    so that one can select the relation form the input file.

    // kozeny-Carman relation
    permeabilityFactor_  =  std::pow(((1-initialPorosity_)/(1-this->porosity_)),2)
            * std::pow((this->porosity_/initialPorosity_),3);

    // Verma-Pruess relation
//  permeabilityFactor_  =  100 * std::pow(((this->porosity_/initialPorosity_)-0.9),2);

    // Modified Fair-Hatch relation with final porosity set to 0.2 and E1=1
//  permeabilityFactor_  =  std::pow((this->porosity_/initialPorosity_),3)
//         * std::pow((std::pow((1 - initialPorosity_),2/3))+(std::pow((0.2 - initialPorosity_),2/3)),2)
//         / std::pow((std::pow((1 -this->porosity_),2/3))+(std::pow((0.2 -this->porosity_),2/3)),2);

    //Timur relation with residual water saturation set to 0.001
//    permeabilityFactor_ =  0.136 * (std::pow(this->porosity_,4.4)) / (2000 * (std::pow(0.001,2)));

    //Timur relation1 with residual water saturation set to 0.001
//    permeabilityFactor_ =  0.136 * (std::pow(this->porosity_,4.4)) / (200000 * (std::pow(0.001,2)));


    //Bern. relation
   // permeabilityFactor_ = std::pow((this->porosity_/initialPorosity_),8);

    //Tixier relation with residual water saturation set to 0.001
    //permeabilityFactor_ = (std::pow((250 * (std::pow(this->porosity_,3)) / 0.001),2)) / initialPermeability_;

    //Coates relation with residual water saturation set to 0.001
    //permeabilityFactor_ = (std::pow((100 * (std::pow(this->porosity_,2)) * (1-0.001) / 0.001,2))) / initialPermeability_ ;


    // energy related quantities not contained in the fluid state
    //asImp_().updateEnergy_(priVars, problem,element, fvGeometry, scvIdx, isOldSol);
    }

   /*!
    * \copydoc ImplicitModel::completeFluidState
    * \param isOldSol Specifies whether this is the previous solution or the current one
    */
  static void completeFluidState(const PrimaryVariables& priVars,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 int scvIdx,
                                 FluidState& fluidState,
                                 bool isOldSol = false)

    {
        Scalar t = Implementation::temperature_(priVars, problem, element,fvGeometry, scvIdx);
        fluidState.setTemperature(t);

        int dofIdxGlobal = problem.model().dofMapper().subIndex(element, scvIdx, dofCodim);
        int phasePresence = problem.model().phasePresence(dofIdxGlobal, isOldSol);

        /////////////
        // set the saturations
        /////////////

        Scalar Sg;
        if (phasePresence == nPhaseOnly)
            Sg = 1.0;
        else if (phasePresence == wPhaseOnly) {
            Sg = 0.0;
        }
        else if (phasePresence == bothPhases) {
            if (formulation == plSg)
                Sg = priVars[switchIdx];
            else if (formulation == pgSl)
                Sg = 1.0 - priVars[switchIdx];
            else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        }
        else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        fluidState.setSaturation(nPhaseIdx, Sg);
        fluidState.setSaturation(wPhaseIdx, 1.0 - Sg);

        /////////////
        // set the pressures of the fluid phases
        /////////////

        // calculate capillary pressure
        const MaterialLawParams &materialParams = problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);
        Scalar pc = MaterialLaw::pc(materialParams, 1 - Sg);

        // extract the pressures
        if (formulation == plSg) {
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx]);
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx] + pc);
        }
        else if (formulation == pgSl) {
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx]);
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx] - pc);
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");

        /////////////
        // calculate the phase compositions
        /////////////

    typename FluidSystem::ParameterCache paramCache;

        // now comes the tricky part: calculate phase composition
        if (phasePresence == bothPhases) {
            // both phases are present, phase composition results from
            // the gas <-> liquid equilibrium. This is
            // the job of the "MiscibleMultiPhaseComposition"
            // constraint solver

            // set the known mole fractions in the fluidState so that they
            // can be used by the Miscible2pNcComposition constraint solver
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(wPhaseIdx, compIdx, priVars[compIdx]);
            }

            Miscible2pNCComposition::solve(fluidState,
                                            paramCache,
                                            wPhaseIdx,  //known phaseIdx
                                            /*setViscosity=*/true,
                                            /*setInternalEnergy=*/false);
        }
        else if (phasePresence == nPhaseOnly){

            Dune::FieldVector<Scalar, numComponents> moleFrac;
            Dune::FieldVector<Scalar, numComponents> fugCoeffL;
            Dune::FieldVector<Scalar, numComponents> fugCoeffG;

            for (int compIdx=0; compIdx<numComponents; ++compIdx)
            {
                fugCoeffL[compIdx] = FluidSystem::fugacityCoefficient(fluidState,
                                        paramCache,
                                        wPhaseIdx,
                                        compIdx);
                fugCoeffG[compIdx] = FluidSystem::fugacityCoefficient(fluidState,
                                        paramCache,
                                        nPhaseIdx,
                                        compIdx);
            }
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                moleFrac[compIdx] = (priVars[compIdx]*fugCoeffL[compIdx]*fluidState.pressure(wPhaseIdx))
                /(fugCoeffG[compIdx]*fluidState.pressure(nPhaseIdx));

            moleFrac[wCompIdx] =  priVars[switchIdx];
            Scalar sumMoleFracNotGas = 0;
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                    sumMoleFracNotGas+=moleFrac[compIdx];
            }
            sumMoleFracNotGas += moleFrac[wCompIdx];
            moleFrac[nCompIdx] = 1 - sumMoleFracNotGas;

//          typedef Dune::FieldMatrix<Scalar, numComponents, numComponents> Matrix;
//          typedef Dune::FieldVector<Scalar, numComponents> Vector;


            // Set fluid state mole fractions
            for (int compIdx=0; compIdx<numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(nPhaseIdx, compIdx, moleFrac[compIdx]);
            }

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase2pNc" constraint solver
            ComputeFromReferencePhase2pNCMin::solve(fluidState,
                                                    paramCache,
                                                    nPhaseIdx,
                                                    nPhaseOnly,
                                                    /*setViscosity=*/true,
                                                    /*setInternalEnergy=*/false);

            }
        else if (phasePresence == wPhaseOnly){

        // only the liquid phase is present, i.e. liquid phase
        // composition is stored explicitly.
        // extract _mass_ fractions in the gas phase
            Dune::FieldVector<Scalar, numComponents> moleFrac;

            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                moleFrac[compIdx] = priVars[compIdx];
            }
            moleFrac[nCompIdx] = priVars[switchIdx];
            Scalar sumMoleFracNotWater = 0;
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                    sumMoleFracNotWater+=moleFrac[compIdx];
            }
            sumMoleFracNotWater += moleFrac[nCompIdx];
            moleFrac[wCompIdx] = 1 -sumMoleFracNotWater;

//             convert mass to mole fractions and set the fluid state
            for (int compIdx=0; compIdx<numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(wPhaseIdx, compIdx, moleFrac[compIdx]);
            }

//             calculate the composition of the remaining phases (as
//             well as the densities of all phases). this is the job
//             of the "ComputeFromReferencePhase2pNc" constraint solver
            ComputeFromReferencePhase2pNCMin::solve(fluidState,
                                                    paramCache,
                                                    wPhaseIdx,
                                                    wPhaseOnly,
                                                    /*setViscosity=*/true,
                                                    /*setInternalEnergy=*/false);
        }
        paramCache.updateAll(fluidState);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            Scalar h = Implementation::enthalpy_(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }
    /*!
     * \brief Returns the volume fraction of the precipitate (solid phase)
     * for the given phaseIdx
     *
     * \param phaseIdx the index of the solid phase
     */
    Scalar precipitateVolumeFraction(int phaseIdx) const
    { return precipitateVolumeFraction_[phaseIdx - numPhases]; }

    /*!
     * \brief Returns the inital porosity of the
     * pure, precipitate-free porous medium
     */
    Scalar initialPorosity() const
    { return initialPorosity_;}

    /*!
     * \brief Returns the inital permeability of the
     * pure, precipitate-free porous medium
     */
    Scalar initialPermeability() const
    { return initialPermeability_;}

    /*!
     * \brief Returns the factor for the reduction of the initial permeability
     * due precipitates in the porous medium
     */
    Scalar permeabilityFactor() const
    { return permeabilityFactor_; }

//    /*!
//     * \brief Returns the mole fraction of a component in the phase
//     *
//     * \param phaseIdx the index of the fluid phase
//     * \param compIdx the index of the component
//     */
//    Scalar moleFraction(int phaseIdx, int compIdx) const
//    {
//       return this->fluidState_.moleFraction(phaseIdx, compIdx);
//    }

    /*!
     * \brief Returns the mole fraction of the salinity in the liquid phase
     */
    Scalar moleFracSalinity() const
    {
        return moleFractionSalinity_;
    }

    /*!
     * \brief Returns the salinity (mass fraction) in the liquid phase
     */
    Scalar salinity() const
    {
        return salinity_;
    }

    /*!
     * \brief Returns the density of the phase for all fluid and solid phases
     *
     * \param phaseIdx the index of the fluid phase
     */
    Scalar density(int phaseIdx) const
    {
        if (phaseIdx < numPhases)
            return this->fluidState_.density(phaseIdx);
        else if (phaseIdx >= numPhases)
            return FluidSystem::precipitateDensity(phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }
    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    {
        if (phaseIdx < numPhases)
            return this->fluidState_.molarDensity(phaseIdx);
        else if (phaseIdx >= numPhases)
            return FluidSystem::precipitateMolarDensity(phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the molality of a component in the phase
     *
     * \param phaseIdx the index of the fluid phase
     * \param compIdx the index of the component
     * molality=\frac{n_{component}}{m_{solvent}}
     * =\frac{n_{component}}{n_{solvent}*M_{solvent}}
     * compIdx of the main component (solvent) in the
     * phase is equal to the phaseIdx
     */
     Scalar molality(int phaseIdx, int compIdx) const // [moles/Kg]
    { return this->fluidState_.moleFraction(phaseIdx, compIdx)
                  /(fluidState_.moleFraction(phaseIdx, phaseIdx)
                  * FluidSystem::molarMass(phaseIdx));}

protected:
    friend class TwoPNCVolumeVariables<TypeTag>;
    static Scalar temperature_(const PrimaryVariables &priVars,
                                const Problem& problem,
                                const Element &element,
                                const FVElementGeometry &fvGeometry,
                                int scvIdx)
    {
        return problem.temperatureAtPos(fvGeometry.subContVol[scvIdx].global);
    }

    template<class ParameterCache>
    static Scalar enthalpy_(const FluidState& fluidState,
                            const ParameterCache& paramCache,
                            int phaseIdx)
    {
        return 0;
    }

   /*!
    * \brief Update all quantities for a given control volume.
    *
    * \param priVars The solution primary variables
    * \param problem The problem
    * \param element The element
    * \param fvGeometry Evaluate function with solution of current or previous time step
    * \param scvIdx The local index of the SCV (sub-control volume)
    * \param isOldSol Evaluate function with solution of current or previous time step
    */
        void updateEnergy_(const PrimaryVariables &priVars,
                           const Problem &problem,
                           const Element &element,
                           const FVElementGeometry &fvGeometry,
                           const int scvIdx,
                           bool isOldSol)
        { };

    Scalar precipitateVolumeFraction_[numSPhases];
//     Scalar saturationIdx_[numSPhases];
    Scalar permeabilityFactor_;
    Scalar initialPorosity_;
    Scalar initialPermeability_;
    Scalar minimumPorosity_;
    Scalar sumPrecipitates_;
    Scalar salinity_;
    Scalar moleFractionSalinity_;
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namespace

#endif
