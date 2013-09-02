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
 * \brief Contains the quantities to calculate the energy flux in the
 *        MpNc box model with kinetic energy transfer enabled.
 */
#ifndef DUMUX_MPNC_ENERGY_FLUX_VARIABLES_KINETIC_HH
#define DUMUX_MPNC_ENERGY_FLUX_VARIABLES_KINETIC_HH

#include <dune/common/fvector.hh>

#include <dumux/common/spline.hh>
#include <dumux/implicit/mpnc/mpncfluxvariables.hh>

namespace Dumux
{

template <class TypeTag>
class MPNCFluxVariablesEnergy<TypeTag, /*enableEnergy=*/true, /*kineticEnergyTransfer=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices)  Indices;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum {dim = GridView::dimension};
    enum {numEnergyEqs             = Indices::numPrimaryEnergyVars};

    typedef Dune::FieldVector<CoordScalar, dim>  DimVector;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*!
     * \brief The constructor
     */
    MPNCFluxVariablesEnergy()
    {}
    /*!
     * \brief update
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param face The SCV (sub-control-volume) face
     * \param fluxVars The flux variables
     * \param elemVolVars The volume variables of the current element
     */
    void update(const Problem & problem,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const SCVFace & face,
                const FluxVariables & fluxVars,
                const ElementVolumeVariables & elemVolVars)
    {
        // calculate temperature gradient using finite element
        // gradients
        DimVector tmp ;

        for(int energyEqIdx=0; energyEqIdx<numEnergyEqs; energyEqIdx++)
            temperatureGradient_[energyEqIdx] = 0.;

        for (int idx = 0;
                idx < face.numFap;
                idx++){
            // FE gradient at vertex idx
            const DimVector & feGrad = face.grad[idx];

            for (int energyEqIdx =0; energyEqIdx < numEnergyEqs; ++energyEqIdx){
                // index for the element volume variables
                int volVarsIdx = face.fapIndices[idx];

                tmp = feGrad;
                tmp   *= elemVolVars[volVarsIdx].temperature(energyEqIdx);
                temperatureGradient_[energyEqIdx] += tmp;
            }
        }
    }

    /*!
     * \brief The total heat flux \f$[J/s]\f$ due to heat conduction
     *        of the rock matrix over the sub-control volume's face.
     *
     * \param energyEqIdx The index of the energy equation
     */
    DimVector temperatureGradient(const unsigned int energyEqIdx) const
    {
        return temperatureGradient_[energyEqIdx];
    }


private:
    DimVector temperatureGradient_[numEnergyEqs];
};

} // end namepace

#endif
