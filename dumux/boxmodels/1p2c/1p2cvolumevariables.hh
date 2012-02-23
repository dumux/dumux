// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2011 by Andreas Lauser                               *
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief Quantities required by the single-phase, two-component box
 *        model defined on a vertex.
 */
#ifndef DUMUX_1P2C_VOLUME_VARIABLES_HH
#define DUMUX_1P2C_VOLUME_VARIABLES_HH

#include <dumux/boxmodels/common/boxvolumevariables.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include "1p2cproperties.hh"

namespace Dumux
{

/*!
 * \ingroup OnePTwoCBoxModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are constant within a
 *        finite volume in the single-phase, two-component model.
 */
template <class TypeTag>
class OnePTwoCVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    typedef typename GET_PROP_TYPE(TypeTag, OnePTwoCIndices) Indices;
    enum {
        phaseIdx = Indices::phaseIdx,
        comp0Idx = Indices::comp0Idx,
        comp1Idx = Indices::comp1Idx,

        pressureIdx = Indices::pressureIdx,
        x1Idx = Indices::x1Idx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dimWorld = GridView::dimensionworld };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar,dimWorld> Vector;

public:
    //! The type returned by the fluidState() method
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars A vector containing the primary variables
     * \param problem The considered problem
     * \param element The considered element of the grid
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param scvIdx  The index of the considered sub-control volume
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

        //calculate all secondary variables from the primary variables and store results in fluidstate
        completeFluidState(priVars, problem, element, elemGeom, scvIdx, fluidState_);

        porosity_ = problem.spatialParameters().porosity(element, elemGeom, scvIdx);
        tortuosity_ = problem.spatialParameters().tortuosity(element, elemGeom, scvIdx);
        dispersivity_ = problem.spatialParameters().dispersivity(element, elemGeom, scvIdx);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState_, phaseIdx);

        diffCoeff_ = FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                             paramCache,
                                                             phaseIdx,
                                                             comp0Idx,
                                                             comp1Idx);

        Valgrind::CheckDefined(porosity_);
        Valgrind::CheckDefined(tortuosity_);
        Valgrind::CheckDefined(dispersivity_);
        Valgrind::CheckDefined(diffCoeff_);

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
        Scalar t = Implementation::temperature_(priVars, problem, element,
                                                elemGeom, scvIdx);
        fluidState.setTemperature(t);

        fluidState.setPressure(phaseIdx, priVars[pressureIdx]);

        Scalar x1 = priVars[x1Idx]; //mole or mass fraction of component 1
        if(!useMoles) //mass-fraction formulation
        {
            // convert mass to mole fractions
            Scalar M0 = FluidSystem::molarMass(comp0Idx);
            Scalar M1 = FluidSystem::molarMass(comp1Idx);
            //meanMolarMass if x1_ is a massfraction
            Scalar meanMolarMass = M0*M1/(M1 + x1*(M0 - M1));

            x1 *= meanMolarMass/M1;
        }
        fluidState.setMoleFraction(phaseIdx, comp0Idx, 1 - x1);
        fluidState.setMoleFraction(phaseIdx, comp1Idx, x1);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, phaseIdx);

        Scalar value;
        value = FluidSystem::density(fluidState, paramCache, phaseIdx);
        fluidState.setDensity(phaseIdx, value);
        value = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
        fluidState.setViscosity(phaseIdx, value);
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     */
    Scalar density() const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     */
    Scalar molarDensity() const
    { return fluidState_.molarDensity(phaseIdx);}

    /*!
     * \brief Return mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
     *
     * \param compIdx The index of the component
     */
    Scalar moleFraction(int compIdx) const
    { return fluidState_.moleFraction(phaseIdx, (compIdx==0)?comp0Idx:comp1Idx); }

    /*!
     * \brief Return mass fraction \f$\mathrm{[kg/kg]}\f$ of a component in the phase.
     *
     * \param compIdx The index of the component
     */
    Scalar massFraction(int compIdx) const
    { return fluidState_.massFraction(phaseIdx, (compIdx==0)?comp0Idx:comp1Idx); }

    /*!
     * \brief Return concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param compIdx The index of the component
     */
    Scalar molarity(int compIdx) const
    { return fluidState_.molarity(phaseIdx, (compIdx==0)?comp0Idx:comp1Idx); }

    /*!
     * \brief Return the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     */
    Scalar pressure() const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ in the fluid.
     */
    Scalar diffCoeff() const
    { return diffCoeff_; }

    /*!
     * \brief Return the tortuosity  \f$\mathrm{[-]}\f$ of the streamlines of the fluid.
     */
    Scalar tortuosity() const
    { return tortuosity_; }

    /*!
     * \brief Returns the dispersivity of the fluid's streamlines.
     */
    const Vector &dispersivity() const
    { return dispersivity_; }

    /*!
     * \brief Return temperature \f$\mathrm{[K]}\f$ inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(phaseIdx); }

    /*!
     * \brief Return the dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of a given phase
     *        within the control volume.
     */
    Scalar viscosity() const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Return the average porosity \f$\mathrm{[-]}\f$ within the control volume.
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

    /*!
     * \brief Called by update() to compute the energy related quantities.
     */
    void updateEnergy_(const PrimaryVariables &sol,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &elemGeom,
                       int vertIdx,
                       bool isOldSol)
    { }

    Scalar porosity_;    //!< Effective porosity within the control volume
    Scalar tortuosity_;
    Vector dispersivity_;
    Scalar diffCoeff_;
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}// end namepace

#endif
