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
#ifndef DUMUX_CC_LOCAL_RESIDUAL_HH
#define DUMUX_CC_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/implicit/localresidual.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup CCModel
 * \ingroup CCLocalResidual
 * \brief Element-wise calculation of the residual for models
 *        based on the fully implicit cell-centered scheme.
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class CCLocalResidual : public ImplicitLocalResidual<TypeTag>
{
    using ParentType = ImplicitLocalResidual<TypeTag>;
    friend class ImplicitLocalResidual<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

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

public:
    // copying the local residual class is not a good idea
    CCLocalResidual(const CCLocalResidual &) = delete;

    CCLocalResidual() : ParentType() {}

protected:

    void evalFluxes_(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const ElementBoundaryTypes& bcTypes,
                     const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // calculate the mass flux over the scv faces
        for (auto&& scvf : scvfs(fvGeometry))
        {
            this->residual_[0] += this->asImp_().computeFlux_(element, fvGeometry, elemVolVars, scvf, bcTypes, elemFluxVarsCache[scvf]);
        }
    }

    PrimaryVariables computeFlux_(const Element &element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& elemVolVars,
                                  const SubControlVolumeFace &scvf,
                                  const ElementBoundaryTypes& bcTypes,
                                  const FluxVariablesCache& fluxVarsCache)
    {
        if (!scvf.boundary() /*TODO: || GET_PROP_VALUE(TypeTag, BoundaryReconstruction)*/)
            return this->asImp_().computeFlux(element, fvGeometry, elemVolVars, scvf, fluxVarsCache);
        else
            return PrimaryVariables(0.0);

    }

    void evalBoundary_(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const ElementBoundaryTypes& elemBcTypes,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
            {
                auto bcTypes = this->problem().boundaryTypes(element, scvf);
                this->residual_[0] += evalBoundaryFluxes_(element, fvGeometry, elemVolVars, scvf, bcTypes, elemFluxVarsCache[scvf]);
            }
        }

        // additionally treat mixed D/N conditions in a strong sense
        if (elemBcTypes.hasDirichlet())
        {
            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                    this->asImp_().evalDirichlet_(element, fvGeometry, elemVolVars, scvf);
            }
        }
    }

    /*!
     * \brief Add all fluxes resulting from Neumann, outflow and pure Dirichlet
     *        boundary conditions to the local residual.
     */
    PrimaryVariables evalBoundaryFluxes_(const Element &element,
                                         const FVElementGeometry& fvGeometry,
                                         const ElementVolumeVariables& elemVolVars,
                                         const SubControlVolumeFace &scvf,
                                         const BoundaryTypes& bcTypes,
                                         const FluxVariablesCache& fluxVarsCache)
    {
        // evaluate the Neumann conditions at the boundary face
        PrimaryVariables flux(0);
        if (bcTypes.hasNeumann() /*TODO: && !GET_PROP_VALUE(TypeTag, BoundaryReconstruction)*/)
            flux += this->asImp_().evalNeumannSegment_(element, fvGeometry, elemVolVars, scvf, bcTypes);

        // TODO: evaluate the outflow conditions at the boundary face
        //if (bcTypes.hasOutflow() /*TODO: && !GET_PROP_VALUE(TypeTag, BoundaryReconstruction)*/)
        //    flux += this->asImp_().evalOutflowSegment_(&intersection, bcTypes);

        // evaluate the pure Dirichlet conditions at the boundary face
        if (bcTypes.hasDirichlet() && !bcTypes.hasNeumann())
            flux += this->asImp_().evalDirichletSegment_(element, fvGeometry, elemVolVars, scvf, bcTypes, fluxVarsCache);

        return flux;
    }


    /*!
     * \brief Evaluate Dirichlet conditions on faces that have mixed
     *        Dirichlet/Neumann boundary conditions.
     */
    void evalDirichlet_(const Element &element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace &scvf)
    {
        BoundaryTypes bcTypes = this->problem().boundaryTypes(element, scvf);

        if (bcTypes.hasDirichlet() && bcTypes.hasNeumann())
            this->asImp_().evalDirichletSegmentMixed_(element, fvGeometry, elemVolVars, scvf, bcTypes);
    }

    /*!
     * \brief Add Neumann boundary conditions for a single scv face
     */
    PrimaryVariables evalNeumannSegment_(const Element& element,
                                         const FVElementGeometry& fvGeometry,
                                         const ElementVolumeVariables& elemVolVars,
                                         const SubControlVolumeFace &scvf,
                                         const BoundaryTypes &bcTypes)
    {
        // temporary vector to store the neumann boundary fluxes
        PrimaryVariables flux(0);

        auto neumannFluxes = this->problem().neumann(element, scvf);

        // multiply neumann fluxes with the area and the extrusion factor
        auto&& scv = fvGeometry.scv(scvf.insideScvIdx());
        neumannFluxes *= scvf.area()*elemVolVars[scv].extrusionFactor();

        // add fluxes to the residual
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            if (bcTypes.isNeumann(eqIdx))
                flux[eqIdx] += neumannFluxes[eqIdx];

        return flux;
    }

    /*!
     * \brief Treat Dirichlet boundary conditions in a weak sense for a single
     *        intersection that only has Dirichlet boundary conditions
     */
    PrimaryVariables evalDirichletSegment_(const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVolumeVariables& elemVolVars,
                                           const SubControlVolumeFace &scvf,
                                           const BoundaryTypes &bcTypes,
                                           const FluxVariablesCache& fluxVarsCache)
    {
        // temporary vector to store the Dirichlet boundary fluxes
        PrimaryVariables flux(0);

        auto dirichletFlux = this->asImp_().computeFlux(element, fvGeometry, elemVolVars, scvf, fluxVarsCache);

        // add fluxes to the residual
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            if (bcTypes.isDirichlet(eqIdx))
                flux[eqIdx] += dirichletFlux[eqIdx];

        return flux;
    }

    /*!
     * \brief Treat Dirichlet boundary conditions in a strong sense for a
     *        single intersection that has mixed D/N boundary conditions
     */
    void evalDirichletSegmentMixed_(const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace &scvf,
                                    const BoundaryTypes &bcTypes)
    {
        // temporary vector to store the Dirichlet values
        PrimaryVariables dirichletValues = this->problem().dirichlet(element, scvf);

        // get the primary variables
        const auto& priVars = elemVolVars[scvf.insideScvIdx()].priVars();

        // set Dirichlet conditions in a strong sense
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (bcTypes.isDirichlet(eqIdx))
            {
                auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                this->residual_[0][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];
            }
        }
    }

    /*!
     * \brief Add outflow boundary conditions for a single intersection
     */
    /*template <class IntersectionIterator>
    void evalOutflowSegment_(const IntersectionIterator &isIt,
                             const BoundaryTypes &bcTypes)
    {
        if (this->element_().geometry().type().isCube() == false)
            DUNE_THROW(Dune::InvalidStateException,
                       "for cell-centered models, outflow BCs only work for cubes.");

        // store pointer to the current FVElementGeometry
        const FVElementGeometry *oldFVGeometryPtr = this->fvElemGeomPtr_;

        // copy the current FVElementGeometry to a local variable
        // and set the pointer to this local variable
        FVElementGeometry fvGeometry = this->fvGeometry_();
        this->fvElemGeomPtr_ = &fvGeometry;

        // get the index of the boundary face
        unsigned bfIdx = isIt->indexInInside();
        unsigned oppositeIdx = bfIdx^1;

        // manipulate the corresponding subcontrolvolume face
        SCVFace& boundaryFace = fvGeometry.boundaryFace[bfIdx];

        // set the second flux approximation index for the boundary face
        for (int nIdx = 0; nIdx < fvGeometry.numNeighbors-1; nIdx++)
        {
            // check whether the two faces are opposite of each other
            if (fvGeometry.subContVolFace[nIdx].fIdx == oppositeIdx)
            {
                boundaryFace.j = nIdx+1;
                break;
            }
        }
        boundaryFace.fapIndices[1] = boundaryFace.j;
        boundaryFace.grad[0] *= -0.5;
        boundaryFace.grad[1] *= -0.5;

        // temporary vector to store the outflow boundary fluxes
        PrimaryVariables values;
        Valgrind::SetUndefined(values);

        this->asImp_().computeFlux(values, bfIdx, true);
        values *= this->curVolVars_(0).extrusionFactor();

        // add fluxes to the residual
        Valgrind::CheckDefined(values);
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (bcTypes.isOutflow(eqIdx))
                this->residual_[0][eqIdx] += values[eqIdx];
        }

        // restore the pointer to the FVElementGeometry
        this->fvElemGeomPtr_ = oldFVGeometryPtr;
    }*/
};

}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH
