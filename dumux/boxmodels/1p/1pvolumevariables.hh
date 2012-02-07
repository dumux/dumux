// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Onur Dogan                                        *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief Quantities required by the one-phase box model defined on a vertex.
 */
#ifndef DUMUX_1P_VOLUME_VARIABLES_HH
#define DUMUX_1P_VOLUME_VARIABLES_HH

#include "1pproperties.hh"
#include <dumux/boxmodels/common/boxvolumevariables.hh>

#include <dumux/material/fluidstates/immisciblefluidstate.hh>

namespace Dumux
{

/*!
 * \ingroup OnePBoxModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the one-phase model.
 */
template <class TypeTag>
class OnePVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, OnePIndices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    //! Type of the fluid state
    typedef Dumux::ImmiscibleFluidState<Scalar, FluidSystem> FluidState;

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
        ParentType::update(priVars, problem, element, elemGeom, scvIdx, isOldSol);

        completeFluidState(priVars, problem, element, elemGeom, scvIdx, fluidState_);
        // porosity
        porosity_ = problem.spatialParameters().porosity(element,
                                                         elemGeom,
                                                         scvIdx);

        // energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(priVars, problem, element, elemGeom, scvIdx, isOldSol);
    };

    /*!
     * \@copydoc BoxModel::completeFluidState
     */
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

        fluidState.setPressure(/*phaseIdx=*/0, priVars[Indices::pressureIdx]);

        // saturation in a single phase is always 1 and thus redundant
        // to set. But since we use the fluid state shared by the
        // immiscible multi-phase models, so we have to set it here...
        fluidState.setSaturation(/*phaseIdx=*/0, 1.0);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, /*phaseIdx=*/0);

        Scalar value = FluidSystem::density(fluidState, paramCache,  /*phaseIdx=*/0);
        fluidState.setDensity(/*phaseIdx=*/0, value);

        value = FluidSystem::viscosity(fluidState, paramCache,  /*phaseIdx=*/0);
        fluidState.setViscosity(/*phaseIdx=*/0, value);
    }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     */
    Scalar pressure() const
    { return fluidState_.pressure(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     */
    Scalar density() const
    { return fluidState_.density(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the dynamic viscosity of the fluid within the
     *        control volume.
     */
    Scalar viscosity() const
    { return fluidState_.viscosity(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the fluid state of the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

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
    Scalar porosity_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
