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
 *        using the non-isothermal Stokes model with a staggered grid.
 *
 */
#ifndef DUMUX_STAGGERED_NAVIERSTOKES_NI_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_NAVIERSTOKES_NI_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/implicit/staggered/localresidual.hh>

#include "properties.hh"

namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(EnableEnergyBalanceStokes);
NEW_PROP_TAG(EnableInertiaTerms);
}

/*!
 * \ingroup StaggeredModel
 * \ingroup StaggeredLocalResidual
 * \brief Element-wise calculation of the residual for models
 * 			based on the staggered grid scheme
 *
 * \todo Please doc me more!
 */

// forward declaration
template<class TypeTag, bool enableComponentTransport, bool enableEnergyBalance>
class StaggeredNavierStokesResidualImpl;

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
        massBalanceIdx = Indices::massBalanceIdx,
        energyBalanceIdx = Indices::energyBalanceIdx
    };

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

//    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief Evaluate the cell center storage terms of the non-isothermal Stokes model.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     * \param scv The sub control volume
     * \param volVars The volume variables
     */
    CellCenterPrimaryVariables computeStorageForCellCenter(const SubControlVolume& scv, const VolumeVariables& volVars)
    {
        CellCenterPrimaryVariables storage(0.0);
//        const Scalar density = useMoles? volVars.molarDensity() : volVars.density();

        // compute storage of mass
        storage = ParentType::computeStorageForCellCenter(scv, volVars);

        // compute the storage of energy
        storage[energyBalanceIdx] = volVars.density() * volVars.internalEnergy();

        return storage;
    }


    /*!
     * \brief Evaluate the cell center source terms of the non-isothermal Stokes model.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     * \param element The element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param globalFaceVars The face variables
     * \param scv The sub control volume
     */
    CellCenterPrimaryVariables computeSourceForCellCenter(const Element &element,
                                                          const FVElementGeometry& fvGeometry,
                                                          const ElementVolumeVariables& elemVolVars,
                                                          const GlobalFaceVars& globalFaceVars,
                                                          const SubControlVolume &scv)
    {
        CellCenterPrimaryVariables source(0.0);

        source = this->problem().sourceAtPos(scv.center())[cellCenterIdx][massBalanceIdx];
        source = this->problem().sourceAtPos(scv.center())[cellCenterIdx][energyBalanceIdx];

        return source;
    }
};

} // end namespace

#endif
