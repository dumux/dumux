// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief Volume averaged quantities required by the Richards model.
 */
#ifndef DUMUX_RICHARDS_VOLUME_VARIABLES_HH
#define DUMUX_RICHARDS_VOLUME_VARIABLES_HH

#include "richardsproperties.hh"
#include "richardsfluidstate.hh"

namespace Dumux
{

/*!
 * \ingroup RichardsModel
 * \brief Volume averaged quantities required by the Richards model.
 *
 * This contains the quantities which are are constant within a finite
 * volume in the Richards model
 */
template <class TypeTag>
class RichardsVolumeVariables : public BoxVolumeVariables<TypeTag>
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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),

        pwIdx = Indices::pwIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GridView::template Codim<0>::Entity Element;

public:
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
        ParentType::update(priVars,
                           problem,
                           element,
                           elemGeom,
                           scvIdx,
                           isOldSol);

        this->asImp_().updateTemperature_(priVars,
                                          element,
                                          elemGeom,
                                          scvIdx,
                                          problem);

        // material law parameters
        Scalar pnRef = problem.referencePressure(element, elemGeom, scvIdx);
        const MaterialLawParams &materialParams =
            problem.spatialParameters().materialLawParams(element, elemGeom, scvIdx);
        fluidState_.update(pnRef, materialParams, priVars, temperature_);

        mobility_[wPhaseIdx] =
            MaterialLaw::krw(materialParams, fluidState_.saturation(wPhaseIdx)) /
            FluidSystem::phaseViscosity(wPhaseIdx,
                                        temperature(),
                                        fluidState_.phasePressure(wPhaseIdx),
                                        fluidState_);
        mobility_[nPhaseIdx] =
            MaterialLaw::krn(materialParams, fluidState_.saturation(wPhaseIdx)) /
            FluidSystem::phaseViscosity(nPhaseIdx,
                                        temperature(),
                                        fluidState_.phasePressure(wPhaseIdx),
                                        fluidState_);
        porosity_ = problem.spatialParameters().porosity(element, elemGeom, scvIdx);
    }

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
    { return fluidState_.phasePressure(phaseIdx); }

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
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the effective capillary pressure \f$\mathrm{[Pa]}\f$ within the
     *        control volume.
     *
     * The capillary pressure is defined as the difference in
     * pressures of the non-wetting and the wetting phase, i.e.
     * \f[ p_c = p_n - p_w \f]
     */
    Scalar capillaryPressure() const
    { return fluidState_.capillaryPressure(); }


protected:
    void updateTemperature_(const PrimaryVariables &priVars,
                            const Element &element,
                            const FVElementGeometry &elemGeom,
                            int scvIdx,
                            const Problem &problem)
    {
        temperature_ = problem.temperature(element, elemGeom, scvIdx);
    }

    FluidState fluidState_;

    Scalar mobility_[numPhases];
    Scalar porosity_;
    Scalar temperature_;
};

}

#endif
