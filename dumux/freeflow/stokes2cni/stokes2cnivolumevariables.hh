// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Klaus Mosthaf                                     *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
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
 *        finite volume in the non-isothermal compositional Stokes model.
 */
#ifndef DUMUX_STOKES2CNI_VOLUME_VARIABLES_HH
#define DUMUX_STOKES2CNI_VOLUME_VARIABLES_HH

#include <dumux/freeflow/stokes2c/stokes2cvolumevariables.hh>
#include "stokes2cniproperties.hh"

namespace Dumux
{

/*!
 * \ingroup BoxStokes2cniModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the non-isothermal two-component Stokes
 *        model.
 */
template <class TypeTag>
class Stokes2cniVolumeVariables : public Stokes2cVolumeVariables<TypeTag>
{
    typedef Stokes2cVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, Stokes2cniIndices) Indices;
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIndex) };
    enum { energyIdx = Indices::energyIdx };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

public:
    /*!
     * \brief Update all additional quantities for a given control volume.
     */
    void update(const PrimaryVariables &primaryVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int vertIdx,
                bool isOldSol)
    {
        // the internal energies and the enthalpies
        heatConductivity_ = 0.0257; //TODO: value (Source: www.engineeringtoolbox.com/air-properties-d_156.html)

        // vertex update data for the mass balance
        ParentType::update(primaryVars,
                           problem,
                           element,
                           elemGeom,
                           vertIdx,
                           isOldSol);
    };

    /*!
     * \brief Returns the total internal energy of a phase in the
     *        sub-control volume.
     */
    Scalar internalEnergy() const
    { return this->fluidState_.internalEnergy(phaseIdx); };

    /*!
     * \brief Returns the total enthalpy of a phase in the sub-control
     *        volume.
     */
    Scalar enthalpy() const
    { return this->fluidState_.enthalpy(phaseIdx); };

    /*!
     * \brief Returns the total heat capacity \f$\mathrm{[J/(K*m^3]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar heatConductivity() const
    { return heatConductivity_; };


protected:
    // this method gets called by the parent class. since this method
    // is protected, we are friends with our parent..
    friend class StokesVolumeVariables<TypeTag>;

    static Scalar temperature_(const PrimaryVariables &priVars,
                            const Problem& problem,
                            const Element &element,
                            const FVElementGeometry &elemGeom,
                            int scvIdx)
    {
        return priVars[energyIdx];
    }

    template<class ParameterCache>
    static Scalar enthalpy_(const FluidState& fluidState,
                            const ParameterCache& paramCache,
                            int phaseIdx)
    {
        return FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
    }

    Scalar heatConductivity_;
};

} // end namepace

#endif
