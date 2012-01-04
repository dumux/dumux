// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \ingroup TwoPNIModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are constant within a
 *        finite volume in the non-isothermal two-phase model.
 */
template <class TypeTag>
class TwoPNIVolumeVariables : public TwoPVolumeVariables<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    enum { temperatureIdx = Indices::temperatureIdx };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    /*!
     * \brief Returns the total internal energy of a phase in the
     *        sub-control volume.
     *
     * \param phaseIdx The phase index
     *
     */
    Scalar internalEnergy(int phaseIdx) const
    { return this->fluidState_.internalEnergy(phaseIdx); };

    /*!
     * \brief Returns the total enthalpy of a phase in the sub-control
     *        volume.
     *
     *  \param phaseIdx The phase index
     */
    Scalar enthalpy(int phaseIdx) const
    { return this->fluidState_.enthalpy(phaseIdx); };

    /*!
     * \brief Returns the total heat capacity \f$\mathrm{[J/K*m^3]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar heatCapacity() const
    { return heatCapacity_; };

protected:
    // this method gets called by the parent class. since this method
    // is protected, we are friends with our parent..
    friend class TwoPVolumeVariables<TypeTag>;

    static Scalar temperature_(const PrimaryVariables &priVars,
                            const Problem& problem,
                            const Element &element,
                            const FVElementGeometry &elemGeom,
                            int scvIdx)
    {
        return priVars[Indices::temperatureIdx];
    }

    template<class ParameterCache>
    static Scalar enthalpy_(const FluidState& fluidState,
                            const ParameterCache& paramCache,
                            int phaseIdx)
    {
        return FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
    }

    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const PrimaryVariables &sol,
                       const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &elemGeom,
                       int scvIdx,
                       bool isOldSol)
    {
        // copmute and set the heat capacity of the solid phase
        heatCapacity_ = problem.spatialParameters().heatCapacity(element, elemGeom, scvIdx);
        Valgrind::CheckDefined(heatCapacity_);
    }

    Scalar heatCapacity_;
};

} // end namepace

#endif
