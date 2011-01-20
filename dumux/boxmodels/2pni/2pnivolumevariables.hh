// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
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
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the non-isothermal two-phase model.
 */
#ifndef DUMUX_2PNI_VOLUME_VARIABLES_HH
#define DUMUX_2PNI_VOLUME_VARIABLES_HH

#include <dumux/boxmodels/2p/2pvolumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPNIBoxModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the non-isothermal two-phase model.
 */
template <class TypeTag>
class TwoPNIVolumeVariables : public TwoPVolumeVariables<TypeTag>
{
    typedef TwoPVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    typedef typename GridView::template Codim<0>::Entity Element;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    enum { temperatureIdx = Indices::temperatureIdx };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars The local primary variable vector
     * \param problem The problem object
     * \param element The current element
     * \param elemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local index of the SCV (sub-control volume)
     * \param isOldSol Evaluate function with solution of current or previous time step
     *
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                bool isOldSol)
    {
        typedef Indices I;

        // vertex update data for the mass balance
        ParentType::update(priVars,
                           problem,
                           element,
                           elemGeom,
                           scvIdx,
                           isOldSol);

        // the internal energies and the enthalpies
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            enthalpy_[phaseIdx] =
                FluidSystem::phaseEnthalpy(phaseIdx,
                                      this->fluidState().temperature(),
                                      this->fluidState().phasePressure(phaseIdx),
                                      this->fluidState());
            internalEnergy_[phaseIdx] =
                FluidSystem::phaseInternalEnergy(phaseIdx,
                                            this->fluidState().temperature(),
                                            this->fluidState().phasePressure(phaseIdx),
                                            this->fluidState());
        }
        Valgrind::CheckDefined(internalEnergy_);
        Valgrind::CheckDefined(enthalpy_);
    }

    /*!
        * \brief Update the temperature for a given control volume.
        *
        * \param priVars The local primary variable vector
        * \param element The current element
        * \param elemGeom The finite-volume geometry in the box scheme
        * \param scvIdx The local index of the SCV (sub-control volume)
        * \param problem The problem object
        *
        */

    // this method gets called by the parent class
    void updateTemperature_(const PrimaryVariables &priVars,
                            const Element &element,
                            const FVElementGeometry &elemGeom,
                            int scvIdx,
                            const Problem &problem)
    {
        // retrieve temperature from primary variables
        this->temperature_ = priVars[temperatureIdx];

        heatCapacity_ =
            problem.spatialParameters().heatCapacity(element, elemGeom, scvIdx);

        Valgrind::CheckDefined(this->temperature_);
        Valgrind::CheckDefined(heatCapacity_);
    }

    /*!
     * \brief Returns the total internal energy of a phase in the
     *        sub-control volume.
     *
     * \param phaseIdx The phase index
     *
     */
    Scalar internalEnergy(int phaseIdx) const
    { return internalEnergy_[phaseIdx]; };

    /*!
     * \brief Returns the total enthalpy of a phase in the sub-control
     *        volume.
     *
     *  \param phaseIdx The phase index
     */
    Scalar enthalpy(int phaseIdx) const
    { return enthalpy_[phaseIdx]; };

    /*!
     * \brief Returns the total heat capacity \f$\mathrm{[J/K*m^3]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar heatCapacity() const
    { return heatCapacity_; };

protected:
    Scalar internalEnergy_[numPhases];
    Scalar enthalpy_[numPhases];
    Scalar heatCapacity_;
};

} // end namepace

#endif
