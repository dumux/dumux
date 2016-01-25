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
 * \brief Contains the energy part of volume variables of the MpNc model.
 */
#ifndef DUMUX_MPNC_ENERGY_VOLUME_VARIABLES_HH
#define DUMUX_MPNC_ENERGY_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/mpnc/implicit/properties.hh>
#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux
{
/*!
 * \brief Contains the energy related quantities which are constant within a
 *        finite volume in a MpNc model.
 *
 * This is the dummy class for the isothermal case. Note that we're
 * only isothermal in the sense that the temperature at a location and
 * a time is specified outside of the model!
 */
template <class TypeTag, bool enableEnergy/*=false*/, int numEnergyEquations /*=don't care*/>
class MPNCVolumeVariablesEnergy
{
    static_assert(not (numEnergyEquations and not enableEnergy),
                  "No kinetic energy transfer may only be enabled "
                  "if energy is enabled in general.");
    static_assert(numEnergyEquations < 1,
                  "No kinetic energy transfer module included, "
                  "but kinetic energy transfer enabled.");

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    enum {dimWorld=GridView::dimensionworld};
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;


public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     *
     * \param fs Container for all the secondary variables concerning the fluids
     * \param paramCache Container for cache parameters
     * \param priVars The primary Variables
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param scvIdx The index of the sub-control volume
     * \param problem The problem
     */
    void updateTemperatures(FluidState &fs,
                            ParameterCache &paramCache,
                            const PrimaryVariables &priVars,
                            const Element &element,
                            const FVElementGeometry &fvGeometry,
                            const unsigned int scvIdx,
                            const Problem &problem) const
    {
        Scalar T = problem.temperatureAtPos(fvGeometry.subContVol[scvIdx].global);
        fs.setTemperature(T);
    }


    /*!
     * \brief Update the enthalpy and the internal energy for a given
     *        control volume.
     * \param fs Container for all the secondary variables concerning the fluids
     * \param paramCache Container for cache parameters
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param scvIdx The index of the sub-control volume
     * \param problem The problem
     *
     * Since we are isothermal, we don't need to do anything!
     */
    void update(FluidState &fs,
                ParameterCache &paramCache,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const unsigned int scvIdx,
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
class MPNCVolumeVariablesEnergy<TypeTag, /*enableEnergy=*/true, /*numEnergyEquations=*/ 1 >
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numPhases        = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { temperatureIdx   = Indices::temperatureIdx };
    enum { numEnergyEqs     = Indices::numPrimaryEnergyVars};
    enum { temperature0Idx = Indices::temperatureIdx };

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState)  FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;

    enum { dimWorld = GridView::dimensionworld};
    typedef Dune::FieldVector<typename GridView::Grid::ctype, dimWorld> GlobalPosition;

public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     *
     * \param fs Container for all the secondary variables concerning the fluids
     * \param paramCache Container for cache parameters
     * \param sol The primary Vaiables
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param scvIdx The index of the sub-control volume
     * \param problem The problem
     */
    void updateTemperatures(FluidState &fs,
                            ParameterCache &paramCache,
                            const PrimaryVariables &sol,
                            const Element &element,
                            const FVElementGeometry &fvGeometry,
                            const unsigned int scvIdx,
                            const Problem &problem) const
    {
        // retrieve temperature from solution vector
        Scalar T = sol[temperatureIdx];
        fs.setTemperature(T);
    }

    /*!
     * \brief Update the enthalpy and the internal energy for a given
     *        control volume.
     * \param fs Container for all the secondary variables concerning the fluids
     * \param paramCache Container for cache parameters
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param scvIdx The index of the sub-control volume
     * \param problem The problem
     */
    void update(FluidState &fs,
                ParameterCache &paramCache,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const unsigned int scvIdx,
                const Problem &problem)
    {
        Valgrind::SetUndefined(*this);

        // heat capacities of the fluids plus the porous medium
        solidHeatCapacity_ =
            problem.spatialParams().solidHeatCapacity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(solidHeatCapacity_);

        solidDensity_ =
            problem.spatialParams().solidDensity(element, fvGeometry, scvIdx);
        Valgrind::CheckDefined(solidDensity_);

        // set the enthalpies
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar h = FluidSystem::enthalpy(fs, paramCache, phaseIdx);
            fs.setEnthalpy(phaseIdx, h);
            thermalConductivity_[phaseIdx] = FluidSystem::thermalConductivity(fs, phaseIdx);
        }
    }

    /*!
     * \brief Returns the total heat capacity [J/(kg K)] of the rock matrix in
     *        the sub-control volume.
     */
    Scalar solidHeatCapacity() const
    { return solidHeatCapacity_; }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ of the fluid phase in
     *        the sub-control volume.
     */
    Scalar thermalConductivity(const unsigned int phaseIdx) const
    { return thermalConductivity_[phaseIdx] ; }

    /*!
     * \brief Returns the total density of the given solid phase [kg / m^3] in
     *        the sub-control volume.
     */
    Scalar solidDensity() const
    { return solidDensity_; }

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    {
        Valgrind::CheckDefined(solidHeatCapacity_);
        Valgrind::CheckDefined(solidDensity_);
    }

protected:
    Scalar solidHeatCapacity_;
    Scalar solidDensity_;
    Scalar thermalConductivity_[numPhases] ;
};

} // end namespace

#endif
