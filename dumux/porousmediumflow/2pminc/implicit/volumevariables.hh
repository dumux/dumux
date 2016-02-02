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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase model.
 */
#ifndef DUMUX_TWOPMINC_VOLUME_VARIABLES_HH
#define DUMUX_TWOPMINC_VOLUME_VARIABLES_HH

//can't be derived from 2pvolumevariables
#include <dumux/implicit/volumevariables.hh>

#include <dune/common/fvector.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPBoxModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase model.
 */
template <class TypeTag>
class TwoPMincVolumeVariables : public ImplicitVolumeVariables<TypeTag>
{
    typedef ImplicitVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pwsn = Indices::pwsn,
        pnsw = Indices::pnsw,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numContinua = GET_PROP_VALUE(TypeTag, NumContinua),
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
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

        completeFluidState(priVars, problem, element, fvGeometry, scvIdx, fluidState_);

        const MaterialLawParams &materialParamsFracture =
            problem.spatialParams().materialLawParamsFracture(element, fvGeometry, scvIdx);
        const MaterialLawParams &materialParamsMatrix =
            problem.spatialParams().materialLawParamsMatrix(element, fvGeometry, scvIdx);

        // relative permeabilities krw/krn for fractures (idx 0) and matrix elements (>= idx 1)
        for (int cIdx = 0; cIdx < numContinua; ++cIdx)
        {
            Scalar krw;
            Scalar krn;
            if (cIdx == 0)
            {
                krw = MaterialLaw::krw(materialParamsFracture,
                        fluidState_[cIdx].saturation(wPhaseIdx));
                krn = MaterialLaw::krn(materialParamsFracture,
                        fluidState_[cIdx].saturation(wPhaseIdx));
            }
            else
            {
                krw = MaterialLaw::krw(materialParamsMatrix,
                        fluidState_[cIdx].saturation(wPhaseIdx));
                krn = MaterialLaw::krn(materialParamsMatrix,
                        fluidState_[cIdx].saturation(wPhaseIdx));
            }

            volFraction_[cIdx] = problem.model().getVolFraction()[cIdx];

            mobility_[wPhaseIdx][cIdx] =
                    krw / fluidState_[cIdx].viscosity(wPhaseIdx);

            mobility_[nPhaseIdx][cIdx] =
                    krn / fluidState_[cIdx].viscosity(nPhaseIdx);
            // porosity
            porosity_[cIdx] = problem.spatialParams().
                    porosity(element,fvGeometry, scvIdx, cIdx);

            intrinsicPermeability_[cIdx] = problem.spatialParams().
                    intrinsicPermeability(element, fvGeometry, scvIdx, cIdx);

            //transferTerm_[cIdx] = problem.spatialParams().
            //        transferTerm(element, fvGeometry, scvIdx, cIdx);

        }

        // energy related quantities not belonging to the fluid state
        asImp_().updateEnergy_(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
    }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     */
    template<typename FluidStateArray>
    static void completeFluidState(const PrimaryVariables &priVars,
                                   const Problem &problem,
                                   const Element &element,
                                   const FVElementGeometry &fvGeometry,
                                   int scvIdx,
                                   FluidStateArray &fluidState)
    {
        Scalar t = Implementation::temperature_(priVars, problem, element,
                                                fvGeometry, scvIdx);
        for (int nC=0; nC<numContinua; nC++)
        {
            fluidState[nC].setTemperature(t);
        }

        const MaterialLawParams &materialParamsFracture =
            problem.spatialParams().materialLawParamsFracture(element, fvGeometry, scvIdx);

        const MaterialLawParams &materialParamsMatrix =
            problem.spatialParams().materialLawParamsMatrix(element, fvGeometry, scvIdx);

        if (int(formulation) == pwsn)
        {
            for (int nC=0; nC<numContinua; nC++)
            {

                Scalar sn = priVars[Indices::sIdxc(nC)];
                fluidState[nC].setSaturation(nPhaseIdx,sn);
                fluidState[nC].setSaturation(wPhaseIdx, 1.0 - sn);

                Scalar pw = priVars[Indices::pIdxc(nC)];
                fluidState[nC].setPressure(wPhaseIdx,pw);

                Scalar pc;
                if (nC == 0){
                    pc = MaterialLaw::pc(materialParamsFracture,
                            1.0-sn);
                }

                else{
                    pc = MaterialLaw::pc(materialParamsMatrix,
                            1.0-sn);
                }
                Scalar pn    = pw + pc;
                fluidState[nC].setPressure(nPhaseIdx,pn);
            }
        }
        else if (int(formulation) == pnsw)
        {
            for (int nC=0; nC<numContinua; nC++)
            {
                Scalar sw = priVars[Indices::sIdxc(nC)];
                fluidState[nC].setSaturation(wPhaseIdx,sw);
                fluidState[nC].setSaturation(nPhaseIdx, 1.0 - sw);

                Scalar pn=priVars[Indices::pIdxc(nC)] ;
                fluidState[nC].setPressure(nPhaseIdx,pn);
                Scalar pc;
                if (nC == 0)
                    pc = MaterialLaw::pc(materialParamsFracture, sw);
                else
                    pc = MaterialLaw::pc(materialParamsMatrix, sw);
                Scalar pw    = pn - pc;
                fluidState[nC].setPressure(wPhaseIdx,pw);
            }
        }

        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typename FluidSystem::ParameterCache paramCache;
        for (int nC=0; nC<numContinua; nC++)
        {
            paramCache.updateAll(fluidState[nC]);
        }


        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int nC=0; nC < numContinua; ++nC)
            {
                // compute and set the viscosity
                Scalar mu = FluidSystem::viscosity(fluidState[nC], paramCache, phaseIdx);
                fluidState[nC].setViscosity(phaseIdx, mu);

                // compute and set the density
                Scalar rho = FluidSystem::density(fluidState[nC], paramCache, phaseIdx);
                fluidState[nC].setDensity(phaseIdx, rho);

                // compute and set the enthalpy
                Scalar h = Implementation::enthalpy_(fluidState[nC], paramCache, phaseIdx);
                fluidState[nC].setEnthalpy(phaseIdx, h);
            }
        }
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     *
     * \param nC The continuum index
     */
    const FluidState &fluidState(int nC) const
    { return fluidState_[nC]; }

    /*!
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     * \param nC The continuum index
     */
    Scalar saturation(int phaseIdx, int nC) const
    { return fluidState_[nC].saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     * \param nC The continuum index
     */
    Scalar density(int phaseIdx, int nC) const
    { return fluidState_[nC].density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     * \param nC The continuum index
     */
    Scalar pressure(int phaseIdx, int nC) const
    { return fluidState_[nC].pressure(phaseIdx); }

    /*!
     * \brief Returns the capillary pressure within the control volume [Pa].
     *
     * \param nC The continuum index
     */
    Scalar capillaryPressure(int nC) const
    { return fluidState_[nC].pressure(nPhaseIdx) - fluidState_[nC].pressure(wPhaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_[/*nC=*/0].temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     * \param nC The continuum index
     */
    Scalar mobility(int phaseIdx, int nC) const
    { return mobility_[phaseIdx][nC]; }

    /*!
     * \brief Returns the average porosity within the control volume.
     *
     * \param nC The continuum index
     */
    Scalar porosity(int nC) const
    { return porosity_[nC]; }

    Scalar intrinsicPermeability(int nC) const
    { return intrinsicPermeability_[nC]; }


    Scalar volumeFraction(int nC) const
    { return volFraction_[nC]; }

protected:
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
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const PrimaryVariables &sol,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &fvGeometry,
                       int vIdx,
                       bool isOldSol)
    { }

    Dune::FieldMatrix<Scalar, numPhases, numContinua> mobility_;
    Dune::FieldVector<Scalar, numContinua> intrinsicPermeability_;

    FluidState fluidState_[numContinua];
    Scalar porosity_[numContinua];
    Scalar volFraction_[numContinua];

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
