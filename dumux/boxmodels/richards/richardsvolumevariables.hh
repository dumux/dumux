// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \brief Volume averaged quantities required by the Richards model.
 */
#ifndef DUMUX_RICHARDS_VOLUME_VARIABLES_HH
#define DUMUX_RICHARDS_VOLUME_VARIABLES_HH

#include "richardsproperties.hh"

#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/boxmodels/common/boxvolumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup RichardsModel
 * \ingroup BoxVolumeVariables
 * \brief Volume averaged quantities required by the Richards model.
 *
 * This contains the quantities which are are constant within a finite
 * volume in the Richards model
 */
template <class TypeTag>
class RichardsVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, RichardsIndices) Indices;
    enum {
        pwIdx = Indices::pwIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    //! The type returned by the fluidState() method
    typedef Dumux::ImmiscibleFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars The primary variables as a vector for the finite
     *                volume.
     * \param problem The physical problem which needs to be solved.
     * \param element The DUNE Codim<0> enitity which intersects
     *                the control volume of the box method
     * \param elemGeom The element's finite volume geometry
     * \param scvIdx The local index of the sub control volume inside the element
     * \param isOldSol Specifies whether the solution is from
     *                 the previous time step or from the current one
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                bool isOldSol)
    {
        assert(!FluidSystem::isLiquid(nPhaseIdx));

        ParentType::update(priVars,
                           problem,
                           element,
                           elemGeom,
                           scvIdx,
                           isOldSol);

        completeFluidState(priVars, problem, element, elemGeom, scvIdx, fluidState_);

        //////////
        // specify the other parameters
        //////////
        const MaterialLawParams &matParams =
            problem.spatialParameters().materialLawParams(element, elemGeom, scvIdx);
        relativePermeabilityWetting_ =
            MaterialLaw::krw(matParams,
                             fluidState_.saturation(wPhaseIdx));

        porosity_ = problem.spatialParameters().porosity(element, elemGeom, scvIdx);

        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(priVars, problem, element, elemGeom, scvIdx, isOldSol);
    }

    /*!
     * \copydoc BoxModel::completeFluidState
     */
    static void completeFluidState(const PrimaryVariables& primaryVariables,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& elementGeometry,
                                   int scvIdx,
                                   FluidState& fluidState)
    {
        // temperature
        Scalar t = Implementation::temperature_(primaryVariables, problem, element,
                                                elementGeometry, scvIdx);
        fluidState.setTemperature(t);

        // pressures
        Scalar pnRef = problem.referencePressure(element, elementGeometry, scvIdx);
        const MaterialLawParams &matParams =
            problem.spatialParameters().materialLawParams(element, elementGeometry, scvIdx);
        Scalar minPc = MaterialLaw::pC(matParams, 1.0);
        fluidState.setPressure(wPhaseIdx, primaryVariables[pwIdx]);
        fluidState.setPressure(nPhaseIdx, std::max(pnRef, primaryVariables[pwIdx] + minPc));

        // saturations
        Scalar Sw = MaterialLaw::Sw(matParams, fluidState.pressure(nPhaseIdx) - fluidState.pressure(wPhaseIdx));
        fluidState.setSaturation(wPhaseIdx, Sw);
        fluidState.setSaturation(nPhaseIdx, 1 - Sw);

        // density and viscosity
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);
        fluidState.setDensity(wPhaseIdx, FluidSystem::density(fluidState, paramCache, wPhaseIdx));
        fluidState.setDensity(nPhaseIdx, 1e-10);

        fluidState.setViscosity(wPhaseIdx, FluidSystem::viscosity(fluidState, paramCache, wPhaseIdx));
        fluidState.setViscosity(nPhaseIdx, 1e-10);
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the average porosity [] within the control volume.
     *
     * The porosity is defined as the ratio of the pore space to the
     * total volume, i.e. \f[ \Phi := \frac{V_{pore}}{V_{pore} + V_{rock}} \f]
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the average absolute saturation [] of a given
     *        fluid phase within the finite volume.
     *
     * The saturation of a fluid phase is defined as the fraction of
     * the pore volume filled by it, i.e.
     * \f[ S_\alpha := \frac{V_\alpha}{V_{pore}} = \phi \frac{V_\alpha}{V} \f]
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the average mass density \f$\mathrm{[kg/m^3]}\f$ of a given
     *        fluid phase within the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar density(int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * For the non-wetting phase (i.e. the gas phase), we assume
     * infinite mobility, which implies that the non-wetting phase
     * pressure is equal to the finite volume's reference pressure
     * defined by the problem.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar pressure(int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns average temperature \f$\mathrm{[K]}\f$ inside the control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

    /*!
     * \brief Returns the effective mobility \f$\mathrm{[1/(Pa*s)]}\f$ of a given phase within
     *        the control volume.
     *
     * The mobility of a fluid phase is defined as the relative
     * permeability of the phase (given by the chosen material law)
     * divided by the dynamic viscosity of the fluid, i.e.
     * \f[ \lambda_\alpha := \frac{k_{r\alpha}}{\mu_\alpha} \f]
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar mobility(int phaseIdx) const
    { return relativePermeability(phaseIdx)/fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns relative permeability [-] of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar relativePermeability(int phaseIdx) const
    {
        if (phaseIdx == wPhaseIdx)
            return relativePermeabilityWetting_;
        return 1;
    }

    /*!
     * \brief Returns the effective capillary pressure \f$\mathrm{[Pa]}\f$ within the
     *        control volume.
     *
     * The capillary pressure is defined as the difference in
     * pressures of the non-wetting and the wetting phase, i.e.
     * \f[ p_c = p_n - p_w \f]
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(nPhaseIdx) - fluidState_.pressure(wPhaseIdx); }

protected:
    static Scalar temperature_(const PrimaryVariables &priVars,
                            const Problem& problem,
                            const Element &element,
                            const FVElementGeometry &elemGeom,
                            int scvIdx)
    {
        return problem.boxTemperature(element, elemGeom, scvIdx);
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
    Scalar relativePermeabilityWetting_;
    Scalar porosity_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
