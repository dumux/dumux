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

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        phaseIdx = Indices::phaseIdx,
        phaseCompIdx = Indices::phaseCompIdx,
        transportCompIdx = Indices::transportCompIdx
    };
    //indices of primary variables
    enum{
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar,dim> DimVector;

public:
    //! The type returned by the fluidState() method
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \copydoc BoxVolumeVariables::update
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx,
                const bool isOldSol)
    {
        ParentType::update(priVars, problem, element, fvGeometry, scvIdx, isOldSol);

        //calculate all secondary variables from the primary variables and store results in fluidstate
        completeFluidState(priVars, problem, element, fvGeometry, scvIdx, fluidState_);

        porosity_ = problem.spatialParams().porosity(element, fvGeometry, scvIdx);
        tortuosity_ = problem.spatialParams().tortuosity(element, fvGeometry, scvIdx);
        dispersivity_ = problem.spatialParams().dispersivity(element, fvGeometry, scvIdx);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState_, phaseIdx);

        diffCoeff_ = FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                             paramCache,
                                                             phaseIdx,
                                                             phaseCompIdx,
                                                             transportCompIdx);

        Valgrind::CheckDefined(porosity_);
        Valgrind::CheckDefined(tortuosity_);
        Valgrind::CheckDefined(dispersivity_);
        Valgrind::CheckDefined(diffCoeff_);

        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
    }

    /*!
     * \copydoc BoxModel::completeFluidState
     */
    static void completeFluidState(const PrimaryVariables& priVars,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const int scvIdx,
                                   FluidState& fluidState)
    {
        Scalar T = Implementation::temperature_(priVars, problem, element,
                                                fvGeometry, scvIdx);
        fluidState.setTemperature(T);

        fluidState.setPressure(phaseIdx, priVars[pressureIdx]);

        Scalar x1 = priVars[massOrMoleFracIdx]; //mole or mass fraction of component 1
        if(!useMoles) //mass-fraction formulation
        {
            // convert mass to mole fractions
            Scalar M0 = FluidSystem::molarMass(phaseCompIdx);
            Scalar M1 = FluidSystem::molarMass(transportCompIdx);
            //meanMolarMass if x1_ is a massfraction
            Scalar meanMolarMass = M0*M1/(M1 + x1*(M0 - M1));

            x1 *= meanMolarMass/M1;
        }
        fluidState.setMoleFraction(phaseIdx, phaseCompIdx, 1 - x1);
        fluidState.setMoleFraction(phaseIdx, transportCompIdx, x1);

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
    { return fluidState_.moleFraction(phaseIdx, (compIdx==0)?phaseCompIdx:transportCompIdx); }

    /*!
     * \brief Return mass fraction \f$\mathrm{[kg/kg]}\f$ of a component in the phase.
     *
     * \param compIdx The index of the component
     */
    Scalar massFraction(int compIdx) const
    { return fluidState_.massFraction(phaseIdx, (compIdx==0)?phaseCompIdx:transportCompIdx); }

    /*!
     * \brief Return concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param compIdx The index of the component
     */
    Scalar molarity(int compIdx) const
    { return fluidState_.molarity(phaseIdx, (compIdx==0)?phaseCompIdx:transportCompIdx); }

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
    const DimVector &dispersivity() const
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
                               const FVElementGeometry &fvGeometry,
                               const int scvIdx)
    {
        return problem.boxTemperature(element, fvGeometry, scvIdx);
    }

    /*!
     * \brief Called by update() to compute the energy related quantities.
     */
    void updateEnergy_(const PrimaryVariables &priVars,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &fvGeometry,
                       const int scvIdx,
                       const bool isOldSol)
    { }

    Scalar porosity_;    //!< Effective porosity within the control volume
    Scalar tortuosity_;
    DimVector dispersivity_;
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
