// $Id$
/*****************************************************************************
 *   Copyright (C) 2008,2009 by Klaus Mosthaf,                               *
 *                              Andreas Lauser,                              *
 *                              Bernd Flemisch                               *
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
 *        finite volume in the two-phase, two-component model.
 */
#ifndef DUMUX_2P2C_VOLUME_VARIABLES_HH
#define DUMUX_2P2C_VOLUME_VARIABLES_HH

#include <dumux/boxmodels/common/boxmodel.hh>
#include <dumux/common/math.hh>

#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

#include "2p2cproperties.hh"
#include "2p2cfluidstate.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, two-component model.
 */
template <class TypeTag>
class TwoPTwoCVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLawParams)) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;
    enum {
        dim = GridView::dimension,

        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        formulation = GET_PROP_VALUE(TypeTag, PTAG(Formulation)),

        lCompIdx = Indices::lCompIdx,
        gCompIdx = Indices::gCompIdx,

        lPhaseIdx = Indices::lPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef TwoPTwoCFluidState<TypeTag> FluidState;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars The primary variables
     * \param problem The problem
     * \param element The element
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

        asImp().updateTemperature_(priVars,
                                   element,
                                   elemGeom,
                                   scvIdx,
                                   problem);

        // capillary pressure parameters
        const MaterialLawParams &materialParams =
            problem.spatialParameters().materialLawParams(element, elemGeom, scvIdx);

        int globalVertIdx = problem.model().dofMapper().map(element, scvIdx, dim);
        int phasePresence = problem.model().phasePresence(globalVertIdx, isOldSol);

        // calculate phase state
        fluidState_.update(priVars, materialParams, temperature(), phasePresence);
        Valgrind::CheckDefined(fluidState_);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // Mobilities
            const Scalar mu =
                FluidSystem::phaseViscosity(phaseIdx,
                                            fluidState().temperature(),
                                            fluidState().phasePressure(lPhaseIdx),
                                            fluidState());
            Scalar kr;
            if (phaseIdx == lPhaseIdx)
                kr = MaterialLaw::krw(materialParams, saturation(lPhaseIdx));
            else // ATTENTION: krn requires the liquid saturation
                // as parameter!
                kr = MaterialLaw::krn(materialParams, saturation(lPhaseIdx));
            mobility_[phaseIdx] = kr / mu;
            Valgrind::CheckDefined(mobility_[phaseIdx]);
            
            // binary diffusion coefficents
            diffCoeff_[phaseIdx] =
                FluidSystem::diffCoeff(phaseIdx,
                                       lCompIdx,
                                       gCompIdx,
                                       fluidState_.temperature(),
                                       fluidState_.phasePressure(phaseIdx),
                                       fluidState_);
            Valgrind::CheckDefined(diffCoeff_[phaseIdx]);
        }

        // porosity
        porosity_ = problem.spatialParameters().porosity(element,
                                                         elemGeom,
                                                         scvIdx);
        Valgrind::CheckDefined(porosity_);
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
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    { return fluidState_.density(phaseIdx) / fluidState_.averageMolarMass(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(int phaseIdx) const
    { return fluidState_.phasePressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
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
     * \brief Returns the binary diffusion coefficients for a phase
     */
    Scalar diffCoeff(int phaseIdx) const
    { return diffCoeff_[phaseIdx]; }


protected:

    void updateTemperature_(const PrimaryVariables &priVars,
                            const Element &element,
                            const FVElementGeometry &elemGeom,
                            int scvIdx,
                            const Problem &problem)
    {
        temperature_ = problem.temperature(element, elemGeom, scvIdx);
    }

    Scalar temperature_;     //!< Temperature within the control volume
    Scalar porosity_;        //!< Effective porosity within the control volume
    Scalar mobility_[numPhases];  //!< Effective mobility within the control volume
    Scalar diffCoeff_[numPhases]; //!< Binary diffusion coefficients for the phases
    FluidState fluidState_;

private:
    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namepace

#endif
