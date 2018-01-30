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
 * \brief Calculates the residual of models based on the box scheme element-wise.
 */
#ifndef DUMUX_STAGGERED_NAVIERSTOKES_NC_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_NAVIERSTOKES_NC_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/implicit/staggered/localresidual.hh>

#include "properties.hh"

namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(EnableComponentTransport);
NEW_PROP_TAG(EnableEnergyBalance);
NEW_PROP_TAG(EnableInertiaTerms);
NEW_PROP_TAG(ReplaceCompEqIdx);
}

/*!
 * \ingroup CCModel
 * \ingroup StaggeredLocalResidual
 * \brief Element-wise calculation of the residual for models
 *        based on the fully implicit cell-centered scheme.
 *
 * \todo Please doc me more!
 */


// // forward declaration
template<class TypeTag, bool enableComponentTransport>
class StaggeredNavierStokesResidualImpl;

// specialization for miscible, isothermal flow
template<class TypeTag>
class StaggeredNavierStokesResidualImpl<TypeTag, true> : public StaggeredNavierStokesResidualImpl<TypeTag, false>
{
    using ParentType = StaggeredNavierStokesResidualImpl<TypeTag, false>;
    friend class StaggeredLocalResidual<TypeTag>;
    friend ParentType;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);


    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    enum {
        conti0EqIdx = Indices::conti0EqIdx,
        phaseIdx = Indices::phaseIdx,

        // The index of the component balance equation
        // that gets replaced with the total mass balance
        replaceCompEqIdx = Indices::replaceCompEqIdx
    };

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);
    using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);

    static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);


     /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the immiscible models.
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     * \note This function should not include the source and sink terms.
     * \note The volVars can be different to allow computing
     *       the implicit euler time derivative here
     */
    CellCenterPrimaryVariables computeStorageForCellCenter(const SubControlVolume& scv,
                                                           const VolumeVariables& volVars) const
    {
        CellCenterPrimaryVariables storage(0.0);

        const Scalar density = useMoles ? volVars.molarDensity() : volVars.density(phaseIdx);

        // compute storage term of all components
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            const int eqIdx = conti0EqIdx + compIdx;

            const Scalar massOrMoleFraction = useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx);
            const Scalar s =  density * massOrMoleFraction;

            if (eqIdx != replaceCompEqIdx)
                storage[eqIdx] += s;
        }
            // in case one balance is substituted by the total mass balance
            if(replaceCompEqIdx < numComponents)
                storage[replaceCompEqIdx] = density;

        EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars);

        return storage;
    }

protected:

    /*!
     * \brief Sets a fixed Dirichlet value for a cell (such as pressure) at the boundary.
     *        This is a provisional alternative to setting the Dirichlet value on the boundary directly.
     *
     * \param insideScv The sub control volume
     * \param elemVolVars The current or previous element volVars
     * \param bcTypes The boundary types
     */
    void setFixedCell_(const SubControlVolume& insideScv,
                       const ElementVolumeVariables& elemVolVars,
                       const BoundaryTypes& bcTypes)
    {
        ParentType::setFixedCell_(insideScv, elemVolVars, bcTypes);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            // get equation index
            const auto eqIdx = conti0EqIdx + compIdx;

            // set a fixed mole fraction for cells
            if(eqIdx != conti0EqIdx && bcTypes.isDirichletCell(eqIdx))
            {
                const auto& insideVolVars = elemVolVars[insideScv];
                const Scalar massOrMoleFraction = useMoles ? insideVolVars.moleFraction(phaseIdx, compIdx) : insideVolVars.massFraction(phaseIdx, compIdx);
                this->ccResidual_[eqIdx] = massOrMoleFraction - this->problem().dirichletAtPos(insideScv.center())[cellCenterIdx][eqIdx];
            }
        }

    }
};
}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH
