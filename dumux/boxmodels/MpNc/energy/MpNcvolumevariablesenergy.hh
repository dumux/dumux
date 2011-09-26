/*****************************************************************************
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
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
 * \brief Contains the energy part of volume variables of the M-phase,
 *        N-component model.
 */
#ifndef DUMUX_MPNC_ENERGY_VOLUME_VARIABLES_HH
#define DUMUX_MPNC_ENERGY_VOLUME_VARIABLES_HH

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

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    //typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCEnergyIndices)) EnergyIndices;


public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     */
    Scalar getTemperature(const PrimaryVariables &sol,
                          const Element &element,
                          const FVElementGeometry &elemGeom,
                          int scvIdx,
                          const Problem &problem,
                          const int temperatureIdx) const
    {
        return problem.boxTemperature(element, elemGeom, scvIdx);
    }


    /*!
     * \brief Update the enthalpy and the internal energy for a given
     *        control volume.
     *
     * Since we are isothermal, we don't need to do anything!
     */
    template <class MutableParams>
    void update(MutableParams &mutParams,
                const PrimaryVariables &sol,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                const Problem &problem)
    {
        temperature_ = problem.boxTemperature(element, elemGeom, scvIdx);

        // set the fluid temperatures
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            mutParams.setTemperature(phaseIdx, temperature_);
    }

    /*!
     * \brief Return the temperature of any given phase [K]
     */
    Scalar temperature(int phaseIdx = 0) const
    { return temperature_; }

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    {
        Valgrind::CheckDefined(temperature_);
    }

protected:
    Scalar temperature_;

};

/*!
 * \brief Contains the energy related quantities which are constant within a
 *        finite volume in the two-phase, N-component model.
 */
template <class TypeTag>
class MPNCVolumeVariablesEnergy<TypeTag, /*enableEnergy=*/true, /*kineticEnergyTransfer=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;

    enum { numPhases        = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents    = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };
    enum { temperatureIdx   = Indices::temperatureIdx };
    enum { numEnergyEqs     = Indices::NumPrimaryEnergyVars};
    enum { temperature0Idx = Indices::temperatureIdx };

public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     */
    Scalar getTemperature(const PrimaryVariables &sol,
                          const Element &element,
                          const FVElementGeometry &elemGeom,
                          int scvIdx,
                          const Problem &problem,
                          const int dummy) const
    {
        // retrieve temperature from solution vector
        return sol[temperatureIdx];
    }

    /*!
     * \brief Update the enthalpy and the internal energy for a given
     *        control volume.
     */
    template <class MutableParams>
    void update(MutableParams &mutParams,
                const PrimaryVariables &sol,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                const Problem &problem)
    {
        Valgrind::SetUndefined(*this);

        // set the fluid temperatures
        temperature_ = sol[temperature0Idx];
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            mutParams.setTemperature(phaseIdx, temperature_);

        heatCapacity_ =
            problem.spatialParameters().heatCapacity(element, elemGeom, scvIdx);
        Valgrind::CheckDefined(heatCapacity_);

        soilDensity_ =
                problem.spatialParameters().soilDensity(element, elemGeom, scvIdx);
        Valgrind::CheckDefined(soilDensity_);

        // set the enthalpies
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar h =
                FluidSystem::computeEnthalpy(mutParams, phaseIdx);
            mutParams.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the total heat capacity [J/(K m^3)] of the rock matrix in
     *        the sub-control volume.
     */
    Scalar heatCapacity() const
    { return heatCapacity_; };

    /*!
     * \brief Returns the temperature in fluid / solid phase(s)
     *        the sub-control volume.
     */
    Scalar temperature(int phaseIdx = 0) const
    { return temperature_; }

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
        Valgrind::CheckDefined(temperature_);
    };

protected:
    Scalar heatCapacity_;
    Scalar soilDensity_;
    Scalar temperature_;
};

} // end namepace

#endif
