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
 *
 * \brief Contains the quantities which are constant within a
 *        finite volume in the compositional n-component Stokes model.
 */
#ifndef DUMUX_STOKESNC_VOLUME_VARIABLES_HH
#define DUMUX_STOKESNC_VOLUME_VARIABLES_HH

#include <dumux/freeflow/stokes/stokesvolumevariables.hh>
#include "stokesncproperties.hh"

namespace Dumux
{

/*!
 * \ingroup BoxStokesncModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the n-component Stokes box model.
 */
template <class TypeTag>
class StokesncVolumeVariables : public StokesVolumeVariables<TypeTag>
{
    typedef StokesVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    //component indices
    enum {  transportCompIdx = Indices::transportCompIdx,
            phaseCompIdx = Indices::phaseCompIdx };
    //number of components
	enum {	numComponents = Indices::numComponents };
    //employed phase index
	enum {	phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx) };
    //primary variable indices
	enum {	massOrMoleFracIdx = Indices::massOrMoleFracIdx };
	//equation indices
    enum {  conti0EqIdx = Indices::conti0EqIdx,
            massBalanceIdx = Indices::massBalanceIdx,
            transportEqIdx = Indices::transportEqIdx };

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \copydoc ImplicitVolumeVariables::update()
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx,
                const bool isOldSol)
    {

		// Model is restricted to 2 components when using mass fractions
		if (!useMoles && numComponents>2)
		{
			DUNE_THROW(Dune::NotImplemented, "This model is restricted to 2 components when using mass fractions!\
			                                  To use mole fractions set property UseMoles true ...");
		}

		// set the mole fractions first
        completeFluidState(priVars, problem, element, fvGeometry, scvIdx, this->fluidState(), isOldSol);

		// update vertex data for the mass and momentum balance
        ParentType::update(priVars,
                           problem,
                           element,
                           fvGeometry,
                           scvIdx,
                           isOldSol);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(this->fluidState());

        for (int compIdx=0; compIdx<numComponents; compIdx++)
        {
			if (phaseCompIdx!=compIdx)
			{
				diffCoeff_[compIdx] = FluidSystem::binaryDiffusionCoefficient(this->fluidState(),
                                                             paramCache,
                                                             phaseIdx,
                                                             compIdx,
                                                             phaseCompIdx);
			}
			else
				diffCoeff_[compIdx] = 0.0;

			Valgrind::CheckDefined(diffCoeff_[compIdx]);
		}
    };

    /*!
     * \copydoc ImplicitModel::completeFluidState()
     * \param isOldSol Specifies whether this is the previous solution or the current one
     */
    static void completeFluidState(const PrimaryVariables& priVars,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const int scvIdx,
                                   FluidState& fluidState,
                                   const bool isOldSol = false)
    {

		if(!useMoles) //mass-fraction formulation
		{
			Scalar massOrMoleFrac[numComponents];
			massOrMoleFrac[transportCompIdx] = priVars[massOrMoleFracIdx];
			massOrMoleFrac[phaseCompIdx] = 1.0 - massOrMoleFrac[transportCompIdx];
			// calculate average molar mass of the gas phase
			Scalar M1 = FluidSystem::molarMass(transportCompIdx);
			Scalar M2 = FluidSystem::molarMass(phaseCompIdx);
			Scalar X2 = massOrMoleFrac[phaseCompIdx];
			Scalar avgMolarMass = M1*M2/(M2 + X2*(M1 - M2));

			// convert mass to mole fractions and set the fluid state
			fluidState.setMoleFraction(phaseIdx, transportCompIdx, massOrMoleFrac[transportCompIdx]*avgMolarMass/M1);
			fluidState.setMoleFraction(phaseIdx, phaseCompIdx, massOrMoleFrac[phaseCompIdx]*avgMolarMass/M2);
		}
		else
		{
			Scalar moleFracPhase, sumMoleFrac(0.0);

			//for components
			for (int compIdx=0; compIdx<numComponents; compIdx++)
			{
				if (conti0EqIdx+compIdx != massBalanceIdx)
				{
					sumMoleFrac += priVars[conti0EqIdx+compIdx];
					fluidState.setMoleFraction(phaseIdx, compIdx, priVars[conti0EqIdx+compIdx]);
				}
			}

			//molefraction for the main component (no primary variable)
			moleFracPhase = 1 - sumMoleFrac;
			fluidState.setMoleFraction(phaseIdx, phaseCompIdx, moleFracPhase);
		}
    }

    /*!
     * \brief Returns the mass fraction of a given component in the
     * 		  given fluid phase within the control volume.
     *
     * \param compIdx The component index
     */
    Scalar massFraction(const int compIdx) const
    { return this->fluidState_.massFraction(phaseIdx, compIdx); }

    /*!
     * \brief Returns the mole fraction of a given component in the
     *        given fluid phase within the control volume.
     *
     * \param compIdx The component index
     */
    Scalar moleFraction(const int compIdx) const
    { return this->fluidState_.moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief Returns the binary (mass) diffusion coefficient
     *
     * \param compIdx The component index
     */
    Scalar diffusionCoeff(int compIdx) const
    { return diffCoeff_[compIdx]; }

protected:
    Scalar diffCoeff_[numComponents]; //!<  Diffusion coefficient
};

} // end namespace

#endif
