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
 *        finite volume in the two-phase model.
 */
#ifndef DUMUX_DECOUPLED_2P_VOLUME_VARIABLES_HH
#define DUMUX_DECOUPLED_2P_VOLUME_VARIABLES_HH

#include "dumux/porousmediumflow/2p/implicit/properties.hh"

#include <dumux/implicit/volumevariables.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase model.
 */
template <class TypeTag>
class DecoupledTwoPVolumeVariables : public ImplicitVolumeVariables<TypeTag>
{
    typedef ImplicitVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
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
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    // export type of fluid state for non-isothermal models
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

        completeFluidState(priVars, problem, element, fvGeometry, scvIdx, fluidState_);

        const auto& materialParams =
            problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);

        mobility_[wPhaseIdx] =
            MaterialLaw::krw(materialParams, fluidState_.saturation(wPhaseIdx))
            / fluidState_.viscosity(wPhaseIdx);

        mobility_[nPhaseIdx] =
            MaterialLaw::krn(materialParams, fluidState_.saturation(wPhaseIdx))
            / fluidState_.viscosity(nPhaseIdx);

        initialPorosity_ = problem.spatialParams().porosity(element,
                                                     fvGeometry,
                                                     scvIdx);
        effPorosity_ = initialPorosity_;

//         if(problem.coupled() == true)
//         {
//             if (isOldSol == true)
//             {
//                 effPorosity_ = problem.getEffPorosityOldTimestep(element,
//                                                         fvGeometry,
//                                                         scvIdx);
//             }
//             else
//             {
//                 effPorosity_ = problem.getEffPorosity(element,
//                                                         fvGeometry,
//                                                         scvIdx);
//             }
//         }

//         Scalar idx = problem.elementMapper().index(element);
//         if(idx == 12873)
//         {
//             Scalar effPorosityOldTimestep = problem.getEffPorosityOldTimestep(element,
//                                                             fvGeometry,
//                                                             scvIdx);
//             Scalar effPorosityNew = problem.getEffPorosity(element,
//                                                             fvGeometry,
//                                                             scvIdx);
//             if(effPorosityOldTimestep > 1e-9)
//             {
//                 std::cout << "effPorosityOldTimestep[" << scvIdx << "] is " << effPorosityOldTimestep << std::endl;
//                 std::cout << "effPorosityNew[" << scvIdx << "] is " << effPorosityNew << std::endl;
//             }
//         }

        // energy related quantities not belonging to the fluid state
        asImp_().updateEnergy_(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
    }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     */
    static void completeFluidState(const PrimaryVariables& priVars,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   int scvIdx,
                                   FluidState& fluidState)
    {
        Scalar t = Implementation::temperature_(priVars, problem, element,
                                                fvGeometry, scvIdx);
        fluidState.setTemperature(t);

        const auto& materialParams =
            problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);

        if (int(formulation) == pwsn) {
            Scalar sn = priVars[saturationIdx];
            fluidState.setSaturation(nPhaseIdx, sn);
            fluidState.setSaturation(wPhaseIdx, 1 - sn);

            Scalar pw = priVars[pressureIdx];
            fluidState.setPressure(wPhaseIdx, pw);
            fluidState.setPressure(nPhaseIdx,
                                   pw + MaterialLaw::pc(materialParams, 1 - sn));
        }
        else if (int(formulation) == pnsw) {
            Scalar sw = priVars[saturationIdx];
            fluidState.setSaturation(wPhaseIdx, sw);
            fluidState.setSaturation(nPhaseIdx, 1 - sw);

            Scalar pn = priVars[pressureIdx];
            fluidState.setPressure(nPhaseIdx, pn);
            fluidState.setPressure(wPhaseIdx,
                                   pn - MaterialLaw::pc(materialParams, sw));
        }

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // compute and set the viscosity
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx, mu);

            // compute and set the density
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);

            // compute and set the enthalpy
            Scalar h = Implementation::enthalpy_(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the phase state for the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume in \f$[kg/m^3]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns the capillary pressure within the control volume
     * in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(nPhaseIdx) - fluidState_.pressure(wPhaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume
     * in \f$[K]\f$.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume in \f$[s*m/kg]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the average porosity within the control volume in \f$[-]\f$.
     */
    Scalar effPorosity() const
    { return effPorosity_; }

    /*!
     * \brief Returns the average porosity within the control volume in \f$[-]\f$.
     */
    Scalar initialPorosity() const
    { return initialPorosity_; }

    Scalar deltaVolumetricStrainCurrentIteration() const
    { return deltaVolumetricStrainCurrentIteration_; }

    Scalar deltaVolumetricStrainOldIteration() const
    { return deltaVolumetricStrainOldIteration_; }


    mutable Scalar effPorosity_;
    mutable Scalar deltaVolumetricStrainCurrentIteration_;
    mutable Scalar deltaVolumetricStrainOldIteration_;

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

    FluidState fluidState_;
//     Scalar effPorosity_;
    Scalar mobility_[numPhases];
    Scalar initialPorosity_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
