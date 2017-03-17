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
 *        using the non-isothermal Stokes box model with a staggered grid.
 *
 */
#ifndef DUMUX_STAGGERED_NAVIERSTOKES_NI_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_NAVIERSTOKES_NI_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh> // TODO ?

#include <dumux/common/valgrind.hh> // TODO ?
#include <dumux/implicit/staggered/localresidual.hh>

#include "properties.hh" // TODO ?

namespace Dumux
{

namespace Properties
{
// forward declaration
//NEW_PROP_TAG(EnableComponentTransport);
NEW_PROP_TAG(EnableEnergyBalance);
NEW_PROP_TAG(EnableInertiaTerms);
//NEW_PROP_TAG(ReplaceCompEqIdx);
}

/*!
 * \ingroup CCModel
 * \ingroup StaggeredLocalResidual
 * \brief Element-wise calculation of the residual for models
 * 			based on the fully implicit cell-centered scheme
 *
 * \todo Please doc me more!
 */

// forward declaration
template<class TypeTag, bool enableComponentTransport, bool enableEnergyBalance>
class StaggeredNavierStokesResidualImpl;

template<class TypeTag>
using StaggeredNavierStokesResidual = StaggeredNavierStokesResidualImpl<TypeTag, GET_PROP_VALUE(TypeTag, EnableComponentTransport),
                                                                 GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;


// specialization for immiscible, non-isothermal flow
template<class TypeTag>
class StaggeredNavierStokesResidualImpl<TypeTag, false, true> : public StaggeredNavierStokesResidualImpl<TypeTag, false, false>
{
    friend class StaggeredLocalResidual<TypeTag>;

    using ParentType = StaggeredNavierStokesResidualImpl<TypeTag, false, false>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    enum {
        energyBalanceIdx = Indices::energyBalanceIdx
    };

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief Evaluate the amount the additional quantities to the Stokes model
     *        (energy equation).
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    CellCenterPrimaryVariables computeStorageForCellCenter(const SubControlVolume& scv, const VolumeVariables& volVars) const
    {
        CellCenterPrimaryVariables storage(0.0);
        const Scalar density = useMoles? volVars.molarDensity() : volVars.density();

        // compute storage of mass
        storage = ParentType::computeStorageForCellCenter(scv, volVars);

        // compute the storage of energy
        storage[energyBalanceIdx] = volVars.density(0) * volVars.internalEnergy(0);

        std::cout << "** Subcontrolvolume " << scv.index() << ": energy storage = "
               << storage[energyBalanceIdx] << std::endl;

        return storage;
    }
};

} // end namespace

#endif
