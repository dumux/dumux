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
 *        finite volume in the three-phase model.
 */
#ifndef DUMUX_3P_VOLUME_VARIABLES_HH
#define DUMUX_3P_VOLUME_VARIABLES_HH

#include <dumux/implicit/common/implicitmodel.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include "3pproperties.hh"

namespace Dumux
{

/*!
 * \ingroup ThreePModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in three-phase model.
 */
template <class TypeTag>
class ThreePVolumeVariables : public ImplicitVolumeVariables<TypeTag>
{
    typedef ImplicitVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        dim = GridView::dimension,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),

        wPhaseIdx = Indices::wPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        swIdx = Indices::swIdx,
        snIdx = Indices::snIdx,
        pressureIdx = Indices::pressureIdx
    };

    static const Scalar R; // universal gas constant

    typedef typename GridView::template Codim<0>::Entity Element;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    //! The type of the object returned by the fluidState() method
    typedef Dumux::ImmiscibleFluidState<Scalar, FluidSystem> FluidState;


    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx,
                bool isOldSol)
    {
        ParentType::update(priVars,
                           problem,
                           element,
                           fvGeometry,
                           scvIdx,
                           isOldSol);

        // capillary pressure parameters
        const MaterialLawParams &materialParams =
            problem.spatialParams().materialLawParams(element, fvGeometry, scvIdx);

        Scalar temp = Implementation::temperature_(priVars, problem, element, fvGeometry, scvIdx);
        fluidState_.setTemperature(temp);

        sw_ = priVars[swIdx];
        sn_ = priVars[snIdx];
        sg_ = 1. - sw_ - sn_;

        Valgrind::CheckDefined(sg_);

        fluidState_.setSaturation(wPhaseIdx, sw_);
        fluidState_.setSaturation(gPhaseIdx, sg_);
        fluidState_.setSaturation(nPhaseIdx, sn_);

        /* now the pressures */
        pg_ = priVars[pressureIdx];

        // calculate capillary pressures
        Scalar pcgw = MaterialLaw::pcgw(materialParams, sw_);
        Scalar pcnw = MaterialLaw::pcnw(materialParams, sw_);
        Scalar pcgn = MaterialLaw::pcgn(materialParams, sw_ + sn_);

        Scalar pcAlpha = MaterialLaw::pcAlpha(materialParams, sn_);
        Scalar pcNW1 = 0.0; // TODO: this should be possible to assign in the problem file

        pn_ = pg_- pcAlpha * pcgn - (1.-pcAlpha)*(pcgw - pcNW1);
        pw_ = pn_ - pcAlpha * pcnw - (1.-pcAlpha)*pcNW1;

        fluidState_.setPressure(wPhaseIdx, pw_);
        fluidState_.setPressure(gPhaseIdx, pg_);
        fluidState_.setPressure(nPhaseIdx, pn_);

        Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
        Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
        Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);

        fluidState_.setDensity(wPhaseIdx, rhoW);
        fluidState_.setDensity(gPhaseIdx, rhoG);
        fluidState_.setDensity(nPhaseIdx, rhoN);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // Mobilities
            const Scalar mu =
                FluidSystem::viscosity(fluidState_,
                                       phaseIdx);
            fluidState_.setViscosity(phaseIdx,mu);

            Scalar kr;
            kr = MaterialLaw::kr(materialParams, phaseIdx,
                                 fluidState_.saturation(wPhaseIdx),
                                 fluidState_.saturation(nPhaseIdx),
                                 fluidState_.saturation(gPhaseIdx));
            mobility_[phaseIdx] = kr / mu;
            Valgrind::CheckDefined(mobility_[phaseIdx]);
        }

        // porosity
        porosity_ = problem.spatialParams().porosity(element,
                                                         fvGeometry,
                                                         scvIdx);
        Valgrind::CheckDefined(porosity_);

        // permeability
        permeability_ = problem.spatialParams().intrinsicPermeability(element,
                                                                          fvGeometry,
                                                                          scvIdx);
        Valgrind::CheckDefined(permeability_);

        // energy related quantities not contained in the fluid state
        typename FluidSystem::ParameterCache paramCache;
        asImp_().updateEnergy_(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // compute and set the enthalpy
            Scalar h = asImp_().enthalpy_(fluidState_, paramCache, phaseIdx);
            fluidState_.setEnthalpy(phaseIdx, h);
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
    Scalar saturation(const int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperatures of the rock matrix and of all fluid phases are
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
    Scalar mobility(const int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \brief Returns the effective capillary pressure within the control volume.
     */
    Scalar capillaryPressure() const
    { return fluidState_.capillaryPressure(); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the permeability within the control volume.
     */
    Scalar permeability() const
    { return permeability_; }

protected:

    static Scalar temperature_(const PrimaryVariables &priVars,
                               const Problem &problem,
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
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const PrimaryVariables &priVars,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &fvGeometry,
                       const int scvIdx,
                       bool isOldSol)
    { }

    Scalar sw_, sg_, sn_, pg_, pw_, pn_;

    Scalar porosity_;        //!< Effective porosity within the control volume
    Scalar permeability_;        //!< Effective porosity within the control volume
    Scalar mobility_[numPhases];  //!< Effective mobility within the control volume
    Scalar bulkDensTimesAdsorpCoeff_; //!< the basis for calculating adsorbed NAPL
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

template <class TypeTag>
const typename ThreePVolumeVariables<TypeTag>::Scalar ThreePVolumeVariables<TypeTag>::R = Constants<typename GET_PROP_TYPE(TypeTag, Scalar)>::R;

} // end namespace

#endif
