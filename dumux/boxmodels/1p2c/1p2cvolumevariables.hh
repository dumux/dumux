// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \ingroup OnePTwoCBoxModel
 * \brief Quantities required by the single-phase, two-component box
 *        model defined on a vertex.
 */
#ifndef DUMUX_1P2C_VOLUME_VARIABLES_HH
#define DUMUX_1P2C_VOLUME_VARIABLES_HH

#include "1p2cfluidstate.hh"

#include <dumux/boxmodels/common/boxvolumevariables.hh>

namespace Dumux
{

/*!
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the single-phase, two-component model.
 */
template <class TypeTag>
class OnePTwoCVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),

        dimWorld = GridView::dimensionworld,

        phaseIdx = Indices::phaseIdx,
        comp0Idx = Indices::comp0Idx,
        comp1Idx = Indices::comp1Idx,

        contiEqIdx = Indices::contiEqIdx,
        transEqIdx = Indices::transEqIdx
    };

    typedef OnePTwoCFluidState<TypeTag> FluidState;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dimWorld> Vector;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars A vector containing the primary variables
     * \param problem The considered problem
     * \param element The considered element of the grid
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param scvIdx  The index of the considered subcontrol volume
     * \param isOldSol Evaluate function with solution of current or previous time step
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                bool isOldSol)
    {
        ParentType::update(priVars, problem, element, elemGeom, scvIdx, isOldSol);
		
		asImp().updateTemperature_(priVars,
                                   element,
                                   elemGeom,
                                   scvIdx,
                                   problem);
                                   
        fluidState_.update(priVars, temperature_);
        
        porosity_ = problem.spatialParameters().porosity(element, elemGeom, scvIdx);
        tortuosity_ = problem.spatialParameters().tortuosity(element, elemGeom, scvIdx);
        dispersivity_ = problem.spatialParameters().dispersivity(element, elemGeom, scvIdx);

        viscosity_ = FluidSystem::phaseViscosity(phaseIdx,
                                                 temperature_,
                                                 pressure(),
                                                 fluidState_);
                                                 
        diffCoeff_ = FluidSystem::diffCoeff(phaseIdx,
                                            comp0Idx,
                                            comp1Idx,
                                            temperature_,
                                            pressure(),
                                            *this);

        Valgrind::CheckDefined(porosity_);
        Valgrind::CheckDefined(viscosity_);
        Valgrind::CheckDefined(tortuosity_);
        Valgrind::CheckDefined(dispersivity_);
        Valgrind::CheckDefined(diffCoeff_);
        Valgrind::CheckDefined(fluidState_);
        Valgrind::CheckDefined(*this);
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns density the of the fluid phase.
     */
    Scalar density() const
    { return fluidState_.density(phaseIdx); }

    Scalar molarDensity() const
    { return fluidState_.molarDensity(phaseIdx);}

    /*!
     * \brief Returns mole fraction of a component in the phase
     *
     * \param compIdx The index of the component
     */
    Scalar moleFrac(int compIdx) const
    { return fluidState_.moleFrac(phaseIdx, (compIdx==0)?comp0Idx:comp1Idx); }

    /*!
     * \brief Returns mass fraction of a component in the phase
     * \param compIdx The index of the component
     */
    Scalar massFrac(int compIdx) const
    { return fluidState_.massFrac(phaseIdx, (compIdx==0)?comp0Idx:comp1Idx); }

    /*!
     * \brief Returns concentration of a component in the phase
     * \param compIdx The index of the component
     */
    Scalar concentration(int compIdx) const
    { return fluidState_.concentration(phaseIdx, (compIdx==0)?comp0Idx:comp1Idx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     */
    Scalar pressure() const
    { return fluidState_.phasePressure(phaseIdx); }

    /*!
     * \brief Returns the binary diffusion coefficient in the fluid
     */
    Scalar diffCoeff() const
    { return diffCoeff_; }

    /*!
     * \brief Returns the tortuosity of the streamlines of the fluid.
     */
    Scalar tortuosity() const
    { return tortuosity_; }

    /*!
     * \brief Returns the dispersivity of the fluid's streamlines.
     */
    const Vector &dispersivity() const
    { return dispersivity_; }

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
     * \brief Returns the dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of a given phase
     *        within the control volume.
     */
    Scalar viscosity() const
    { return viscosity_; }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

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
    Scalar porosity_;	//!< Effective porosity within the control volume
    Scalar viscosity_;
    Scalar tortuosity_;
    Vector dispersivity_;
    Scalar diffCoeff_;
    FluidState fluidState_;
    
private:
    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }
};

}// end namepace

#endif
