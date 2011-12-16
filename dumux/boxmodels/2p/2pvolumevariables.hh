/*****************************************************************************
 *   Copyright (C) 2008-2011 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
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
#ifndef DUMUX_2P_VOLUME_VARIABLES_HH
#define DUMUX_2P_VOLUME_VARIABLES_HH

#include "2pproperties.hh"

#include <dumux/boxmodels/common/boxvolumevariables.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPBoxModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase model.
 */
template <class TypeTag>
class TwoPVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLawParams)) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),

        formulation = GET_PROP_VALUE(TypeTag, PTAG(Formulation)),

        pwSn = Indices::pwSn,
        pnSw = Indices::pnSw,

        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GridView::template Codim<0>::Entity Element;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars The local primary variable vector
     * \param problem The problem object
     * \param element The current element
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local index of the SCV (sub-control volume)
     * \param isOldSol Evaluate function with solution of current or previous time step
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                bool isOldSol)
    {
        ParentType::update(priVars,
                           problem,
                           element,
                           elemGeom,
                           scvIdx,
                           isOldSol);

        completeFluidState(priVars, problem, element, elemGeom, scvIdx, fluidState_);

        const MaterialLawParams &materialParams =
            problem.spatialParameters().materialLawParams(element, elemGeom, scvIdx);

        mobility_[wPhaseIdx] =
            MaterialLaw::krw(materialParams, fluidState_.saturation(wPhaseIdx))
            / fluidState_.viscosity(wPhaseIdx);

        mobility_[nPhaseIdx] =
            MaterialLaw::krn(materialParams, fluidState_.saturation(wPhaseIdx))
            / fluidState_.viscosity(nPhaseIdx);

        // porosity
        porosity_ = problem.spatialParameters().porosity(element,
                                                         elemGeom,
                                                         scvIdx);

        // energy related quantities not belonging to the fluid state
        asImp_().updateEnergy_(priVars, problem, element, elemGeom, scvIdx, isOldSol);
    }

    static void completeFluidState(const PrimaryVariables& priVars,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& elemGeom,
                                   int scvIdx,
                                   FluidState& fluidState)
    {
        Scalar t = Implementation::temperature_(priVars, problem, element,
                                                elemGeom, scvIdx);
        fluidState.setTemperature(t);

        // material law parameters
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
        const typename MaterialLaw::Params &materialParams =
            problem.spatialParameters().materialLawParams(element, elemGeom, scvIdx);


        if (int(formulation) == pwSn) {
            Scalar Sn = priVars[saturationIdx];
            fluidState.setSaturation(nPhaseIdx, Sn);
            fluidState.setSaturation(wPhaseIdx, 1 - Sn);

            Scalar pW = priVars[pressureIdx];
            fluidState.setPressure(wPhaseIdx, pW);
            fluidState.setPressure(nPhaseIdx,
                                   pW + MaterialLaw::pC(materialParams, 1 - Sn));
        }
        else if (int(formulation) == pnSw) {
            Scalar Sw = priVars[saturationIdx];
            fluidState.setSaturation(wPhaseIdx, Sw);
            fluidState.setSaturation(nPhaseIdx, 1 - Sw);

            Scalar pN = priVars[pressureIdx];
            fluidState.setPressure(nPhaseIdx, pN);
            fluidState.setPressure(wPhaseIdx,
                                   pN - MaterialLaw::pC(materialParams, Sw));
        }

        typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
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
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns the capillary pressure within the control volume [Pa].
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(nPhaseIdx) - fluidState_.pressure(wPhaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

protected:
    static Scalar temperature_(const PrimaryVariables &priVars,
                            const Problem& problem,
                            const Element &element,
                            const FVElementGeometry &elemGeom,
                            int scvIdx)
    {
        return problem.boxTemperature(element, elemGeom, scvIdx);
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
                       const FVElementGeometry &elemGeom,
                       int vertIdx,
                       bool isOldSol)
    { }

    FluidState fluidState_;
    Scalar porosity_;
    Scalar mobility_[numPhases];

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
