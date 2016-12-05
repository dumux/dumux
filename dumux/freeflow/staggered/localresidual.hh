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
#ifndef DUMUX_STAGGERED_NAVIERSTOKES_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_NAVIERSTOKES_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/implicit/staggered/localresidual.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup CCModel
 * \ingroup StaggeredLocalResidual
 * \brief Element-wise calculation of the residual for models
 *        based on the fully implicit cell-centered scheme.
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class StaggeredNavierStokesResidual : public Dumux::StaggeredLocalResidual<TypeTag>
{
    using ParentType = StaggeredLocalResidual<TypeTag>;
    friend class StaggeredLocalResidual<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using CellCenterSolutionVector = typename GET_PROP_TYPE(TypeTag, CellCenterSolutionVector);
    using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);


    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    enum {
         // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);

public:
    // copying the local residual class is not a good idea
    StaggeredNavierStokesResidual(const StaggeredNavierStokesResidual &) = delete;

    StaggeredNavierStokesResidual() = default;


    static CellCenterPrimaryVariables computeFluxForCellCenter(const Element &element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& elemVolVars,
                                  const GlobalFaceVars& globalFaceVars,
                                  const SubControlVolumeFace &scvf,
                                  const FluxVariablesCache& fluxVarsCache)
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const Scalar velocity = globalFaceVars.faceVars(scvf.dofIndexSelf()).velocity();

        CellCenterPrimaryVariables flux(0.0);

        if(scvf.unitOuterNormal()[scvf.directionIndex()] > 0.0) // positive coordinate direction
        {
            if(velocity > 0.0)
                flux[0] = insideVolVars.density(0) * velocity;
            else
                flux[0] = outsideVolVars.density(0) * velocity;
        }
        else // negative coordinate direction
        {
            if(velocity > 0.0)
                flux[0] = outsideVolVars.density(0) * velocity;
            else
                flux[0] = insideVolVars.density(0) * velocity;
        }
        return flux;
    }

    static CellCenterPrimaryVariables computeSourceForCellCenter(const Element &element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& elemVolVars,
                                     const GlobalFaceVars& globalFaceVars,
                                     const SubControlVolume &scv)
    {
        return CellCenterPrimaryVariables(0.0);
    }


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
    static CellCenterPrimaryVariables computeStorageForCellCenter(const SubControlVolume& scv,
                                    const VolumeVariables& volVars)
    {
        CellCenterPrimaryVariables storage;
        storage[0] = volVars.density(0);
        return storage;
    }

     /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the immiscible models.
     * \param scvf The sub control volume
     * \param volVars The current or previous volVars
     * \note This function should not include the source and sink terms.
     * \note The volVars can be different to allow computing
     *       the implicit euler time derivative here
     */
    static FacePrimaryVariables computeStorageForFace(const SubControlVolumeFace& scvf,
                                    const VolumeVariables& volVars,
                                    const GlobalFaceVars& globalFaceVars)
    {
        FacePrimaryVariables storage;
        const Scalar velocity = globalFaceVars.faceVars(scvf.dofIndexSelf()).velocity();
        storage[0] = volVars.density(0) * velocity;
        return storage;
    }

    static FacePrimaryVariables computeSourceForFace(const SubControlVolumeFace& scvf,
                                    const ElementVolumeVariables& elemVolVars,
                                    const GlobalFaceVars& globalFaceVars)
    {
        return FacePrimaryVariables(0.0);
    }

    static FacePrimaryVariables computeFluxForFace(const SubControlVolumeFace& scvf,
                                    const ElementVolumeVariables& elemVolVars,
                                    const GlobalFaceVars& globalFaceVars)
    {
        return FacePrimaryVariables(0.0);
    }



};

}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH
