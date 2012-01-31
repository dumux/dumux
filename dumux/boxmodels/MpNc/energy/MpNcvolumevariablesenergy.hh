// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
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
 * \brief Contains the energy part of volume variables of the M-phase,
 *        N-component model.
 */
#ifndef DUMUX_MPNC_ENERGY_VOLUME_VARIABLES_HH
#define DUMUX_MPNC_ENERGY_VOLUME_VARIABLES_HH

#include <dumux/boxmodels/MpNc/MpNcproperties.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

namespace Dumux
{
/*!
 * \brief Contains the energy related quantities which are constant within a
 *        finite volume in the two-phase, N-component model.
 *
 * This is the dummy class for the isothermal case. Note that we're
 * only isothermal in the sense that the temperature at a location and
 * a time is specified outside of the model!
 */
template <class TypeTag, bool enableEnergy/*=false*/, bool kineticEnergyTransfer /*=don't care*/>
class MPNCVolumeVariablesEnergy
{
    static_assert(!(kineticEnergyTransfer && !enableEnergy),
                  "No kinetic energy transfer may only be enabled "
                  "if energy is enabled in general.");
    static_assert(!kineticEnergyTransfer,
                  "No kinetic energy transfer module included, "
                  "but kinetic energy transfer enabled.");

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    //typedef typename GET_PROP_TYPE(TypeTag, MPNCEnergyIndices) EnergyIndices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;

public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     */
    void updateTemperatures(FluidState &fs,
                            ParameterCache &paramCache,
                            const PrimaryVariables &sol,
                            const Element &element,
                            const FVElementGeometry &elemGeom,
                            int scvIdx,
                            const Problem &problem) const
    {
        Scalar T = problem.boxTemperature(element, elemGeom, scvIdx);
        fs.setTemperature(T);
    }


    /*!
     * \brief Update the enthalpy and the internal energy for a given
     *        control volume.
     *
     * Since we are isothermal, we don't need to do anything!
     */
    void update(FluidState &fs,
                ParameterCache &paramCache,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                const Problem &problem)
    {
    }

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    {
    }
};

/*!
 * \brief Contains the energy related quantities which are constant within a
 *        finite volume in the two-phase, N-component model.
 */
template <class TypeTag>
class MPNCVolumeVariablesEnergy<TypeTag, /*enableEnergy=*/true, /*kineticEnergyTransfer=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, MPNCIndices) Indices;

    enum { numPhases        = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents    = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { temperatureIdx   = Indices::temperatureIdx };
    enum { numEnergyEqs     = Indices::NumPrimaryEnergyVars};
    enum { temperature0Idx = Indices::temperatureIdx };

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> FluidState;

public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     */
    void updateTemperatures(FluidState &fs,
                            ParameterCache &paramCache,
                            const PrimaryVariables &sol,
                            const Element &element,
                            const FVElementGeometry &elemGeom,
                            int scvIdx,
                            const Problem &problem) const
    {
        // retrieve temperature from solution vector
        Scalar T = sol[temperatureIdx];
        fs.setTemperature(T);
    }

    /*!
     * \brief Update the enthalpy and the internal energy for a given
     *        control volume.
     */
    void update(FluidState &fs,
                ParameterCache &paramCache,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                const Problem &problem)
    {
        Valgrind::SetUndefined(*this);

        // heat capacities of the fluids plus the porous medium
        heatCapacity_ =
            problem.spatialParameters().heatCapacity(element, elemGeom, scvIdx);
        Valgrind::CheckDefined(heatCapacity_);

        soilDensity_ =
            problem.spatialParameters().soilDensity(element, elemGeom, scvIdx);
        Valgrind::CheckDefined(soilDensity_);

        // set the enthalpies
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar h = FluidSystem::enthalpy(fs, paramCache, phaseIdx);
            fs.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the total heat capacity [J/(K m^3)] of the rock matrix in
     *        the sub-control volume.
     */
    Scalar heatCapacity() const
    { return heatCapacity_; };

    /*!
     * \brief Returns the total density of the given soil [kg / m^3] in
     *        the sub-control volume.
     */
    Scalar soilDensity() const
    { return soilDensity_; };

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    {
        Valgrind::CheckDefined(heatCapacity_);
        Valgrind::CheckDefined(soilDensity_);
    };

protected:
    Scalar heatCapacity_;
    Scalar soilDensity_;
};

} // end namepace

#endif
