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
 *        finite volume in the two-phase, one-component model.
 *
 * \note The 2p1c model requires the use of the non-isothermal extension found in dumux/implicit/nonisothermal
 */
#ifndef DUMUX_2P1C_VOLUME_VARIABLES_HH
#define DUMUX_2P1C_VOLUME_VARIABLES_HH

#include <dumux/implicit/volumevariables.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPOneCModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */
template <class TypeTag>
class TwoPOneCVolumeVariables : public ImplicitVolumeVariables<TypeTag>
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
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        dim = GridView::dimension,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        wPhaseIdx = Indices::wPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        switch1Idx = Indices::switch1Idx,
        pressureIdx = Indices::pressureIdx
    };

    // present phases
    enum {
        twoPhases = Indices::twoPhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        gPhaseOnly  = Indices::gPhaseOnly,
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    //! The type of the object returned by the fluidState() method
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;

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

        int globalIdx = problem.model().dofMapper().subIndex(element, scvIdx, dofCodim);
        int phasePresence = problem.model().phasePresence(globalIdx, isOldSol);

        // get saturations
        Scalar sw(0.0);
        Scalar sg(0.0);
        if (phasePresence == twoPhases)
        {
            sw = priVars[switch1Idx];
            sg = 1.0 - sw;
        }
        else if (phasePresence == wPhaseOnly)
        {
            sw = 1.0;
            sg = 0.0;
        }
        else if (phasePresence == gPhaseOnly)
        {
            sw = 0.0;
            sg = 1.0;
        }
        else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        Valgrind::CheckDefined(sg);

        fluidState_.setSaturation(wPhaseIdx, sw);
        fluidState_.setSaturation(gPhaseIdx, sg);

        // get gas phase pressure
        const Scalar pg = priVars[pressureIdx];

        // calculate capillary pressure
        const Scalar pc = MaterialLaw::pc(materialParams, sw);

        // set wetting phase pressure
        const Scalar pw = pg - pc;

        //set pressures
        fluidState_.setPressure(wPhaseIdx, pw);
        fluidState_.setPressure(gPhaseIdx, pg);

        // get temperature
        Scalar temperature;
        if (phasePresence == wPhaseOnly || phasePresence == gPhaseOnly)
            temperature = priVars[switch1Idx];
        else if (phasePresence == twoPhases)
            temperature = FluidSystem::vaporTemperature(fluidState_, wPhaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

        Valgrind::CheckDefined(temperature);

        fluidState_.setTemperature(temperature);

        // set the densities
        const Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
        const Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);

        fluidState_.setDensity(wPhaseIdx, rhoW);
        fluidState_.setDensity(gPhaseIdx, rhoG);

        //get the viscosity and mobility
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // Mobilities
            const Scalar mu =
                FluidSystem::viscosity(fluidState_,
                                       phaseIdx);
            fluidState_.setViscosity(phaseIdx,mu);

            Scalar kr;
            if (phaseIdx == wPhaseIdx)
                kr = MaterialLaw::krw(materialParams, saturation(wPhaseIdx));
            else // ATTENTION: krn requires the wetting phase saturation
                // as parameter!
                kr = MaterialLaw::krn(materialParams, saturation(wPhaseIdx));

            mobility_[phaseIdx] = kr / mu;
            Valgrind::CheckDefined(mobility_[phaseIdx]);
        }

        // porosity
        porosity_ = problem.spatialParams().porosity(element,
                                                         fvGeometry,
                                                         scvIdx);
        Valgrind::CheckDefined(porosity_);

        // the enthalpies (internal energies are directly calculated in the fluidstate
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            Scalar h = FluidSystem::enthalpy(fluidState_, phaseIdx);
            fluidState_.setEnthalpy(phaseIdx, h);
        }

        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
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
     * \brief Returns the molar density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(const int phaseIdx) const
    { return fluidState_.density(phaseIdx) / fluidState_.averageMolarMass(phaseIdx); }

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
     * \brief Returns the vapor temperature (T_{vap}(p_g) of the fluid within the control volume.
     */
    Scalar vaporTemperature() const
    { return FluidSystem::vaporTemperature(fluidState_, wPhaseIdx);}

protected:
    FluidState fluidState_;

private:

    Scalar porosity_;
    Scalar mobility_[numPhases];

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namespace

#endif
