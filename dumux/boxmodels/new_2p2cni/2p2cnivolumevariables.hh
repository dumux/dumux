// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis                                  *
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
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
 * \brief Contains the quantities which are constant within a
 *        finite volume in the non-isothermal two-phase, two-component
 *        model.
 */
#ifndef DUMUX_2P2CNI_VOLUME_VARIABLES_HH
#define DUMUX_2P2CNI_VOLUME_VARIABLES_HH

#include <dumux/boxmodels/new_2p2c/2p2cvolumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCNIModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the non-isothermal two-phase, two-component
 *        model.
 */
template <class TypeTag>
class TwoPTwoCNIVolumeVariables : public TwoPTwoCVolumeVariables<TypeTag>
{
    //! \cond 0
    typedef TwoPTwoCVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };
    enum { temperatureIdx = Indices::temperatureIdx };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    //! \endcond

public:
    /*!
     * \brief Update the temperature of the sub-control volume.
     */
    Scalar getTemperature(const PrimaryVariables &sol,
                          const Element &element,
                          const FVElementGeometry &elemGeom,
                          int scvIdx,
                          const Problem &problem) const
    {
        // retrieve temperature from solution vector
        return sol[temperatureIdx];
    }
    
    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param sol The solution primary variables
     * \param problem The problem
     * \param element The element
     * \param elemGeom Evaluate function with solution of current or previous time step
     * \param vertIdx The local index of the SCV (sub-control volume)
     * \param isOldSol Evaluate function with solution of current or previous time step
     */
    template <class MutableParams>
    void updateEnergy(MutableParams &mutParams,
                      const PrimaryVariables &sol,
                      const Element &element,
                      const FVElementGeometry &elemGeom,
                      int scvIdx,
                      const Problem &problem)
    {
        heatCapacity_ =
            problem.spatialParameters().heatCapacity(element, elemGeom, scvIdx);
        Valgrind::CheckDefined(heatCapacity_);

        // the internal energies and the enthalpies
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar h =
                FluidSystem::computeEnthalpy(mutParams, phaseIdx);
            mutParams.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the total heat capacity \f$\mathrm{[J/(K*m^3]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar heatCapacity() const
    { return heatCapacity_; };

protected:
    Scalar heatCapacity_;
};

} // end namepace

#endif
