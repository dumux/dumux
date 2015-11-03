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
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase n-component mineralisation box model.
 */

#ifndef DUMUX_2PNCMIN_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_2PNCMIN_LOCAL_RESIDUAL_BASE_HH

#include "2pncminproperties.hh"
#include <dumux/implicit/2pnc/2pnclocalresidual.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPNCMinModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase n-component mineralization fully implicit box model.
 *
 * This class is used to fill the gaps in ImplicitLocalResidual for the two-phase n-component flow.
 */
template<class TypeTag>
class TwoPNCMinLocalResidual: public TwoPNCLocalResidual<TypeTag>
{
protected:
    typedef TwoPNCLocalResidual<TypeTag> ParentType;
    typedef TwoPNCMinLocalResidual<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, ElementSolutionVector) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;


    enum
    {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numSPhases = GET_PROP_VALUE(TypeTag, NumSPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx),

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        conti0EqIdx = Indices::conti0EqIdx,

        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases,

        plSg = TwoPNCFormulation::plSg,
        pgSl = TwoPNCFormulation::pgSl,
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    TwoPNCMinLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        this->massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    };

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume).
     * In contrast to the 2pnc model, here, the storage of solid phases is included too.
     *
     *  \param storage the mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
  void computeStorage(PrimaryVariables &storage, int scvIdx, bool usePrevSol) const
  {
      //call parenttype function
      ParentType::computeStorage(storage, scvIdx, usePrevSol);

      const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_()
      : this->curVolVars_();
    const VolumeVariables &volVars = elemVolVars[scvIdx];

    // Compute storage term of all solid (precipitated) phases (excluding the non-reactive matrix)
    for (int phaseIdx = numPhases; phaseIdx < numPhases + numSPhases; ++phaseIdx)
    {
      int eqIdx = conti0EqIdx + numComponents-numPhases + phaseIdx;
      storage[eqIdx] += volVars.precipitateVolumeFraction(phaseIdx)*volVars.molarDensity(phaseIdx);
    }

      Valgrind::CheckDefined(storage);
  }
};
} // end namespace

#endif
