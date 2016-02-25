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
#ifndef DUMUX_CC_TPFA_LOCAL_RESIDUAL_HH
#define DUMUX_CC_TPFA_LOCAL_RESIDUAL_HH

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
class CCTpfaLocalResidual : public ImplicitLocalResidual<TypeTag>
{
    typedef ImplicitLocalResidual<TypeTag> ParentType;
    friend class ImplicitLocalResidual<TypeTag>;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq)
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    // copying the local residual class is not a good idea
    CCTpfaLocalResidual(const CCTpfaLocalResidual &);

public:
    CCTpfaLocalResidual() : ParentType()
    { }

protected:

    void computeFlux_(PrimaryVariables &flux, const SubControlVolumeFace &scvFace)
    {
        if (!scvFace.boundary() /*TODO: || GET_PROP_VALUE(TypeTag, BoundaryReconstruction)*/)
        {
            return this->asImp_().computeFlux(scvFace);
        }
        else
            return this->asImp_().evalBoundary_(scvFace);
    }

    PrimaryVariables evalBoundary_(const SubControlVolumeFace &scvFace)
    {
        PrimaryVariables flux = evalBoundaryFluxes_(scvFace);

        // additionally treat mixed D/N conditions in a strong sense
        if (bcTypes_().hasDirichlet())
            this->asImp_().evalDirichlet_(flux, scvFace);

        return flux;
    }

    /*!
     * \brief Add all fluxes resulting from Neumann, outflow and pure Dirichlet
     *        boundary conditions to the local residual.
     */
    PrimaryVariables evalBoundaryFluxes_(const SubControlVolumeFace &scvFace)
    {
        BoundaryTypes bcTypes = this->problem_().boundaryTypes(this->element_(), scvFace);

        // evaluate the Neumann conditions at the boundary face
        PrimaryVariables flux(0);
        if (bcTypes.hasNeumann())
            flux += this->asImp_().evalNeumannSegment_(scvFace, bcTypes);

        // TODO: evaluate the outflow conditions at the boundary face
        //if (bcTypes.hasOutflow())
        //    flux += this->asImp_().evalOutflowSegment_(&intersection, bcTypes);

        // evaluate the pure Dirichlet conditions at the boundary face
        if (bcTypes.hasDirichlet() && !bcTypes.hasNeumann())
            flux += this->asImp_().evalDirichletSegment_(scvFace, bcTypes);

        return flux;
    }

    /*!
     * \brief Evaluate Dirichlet conditions on faces that have mixed
     *        Dirichlet/Neumann boundary conditions.
     */
    void evalDirichlet_(PrimaryVariables &flux, const SubControlVolumeFace &scvFace)
    {
        BoundaryTypes bcTypes = this->problem_().boundaryTypes(this->element_(), scvFace);

        if (bcTypes.hasDirichlet() && bcTypes.hasNeumann())
            this->asImp_().evalDirichletSegmentMixed_(flux, scvFace);
    }

    /*!
     * \brief Add Neumann boundary conditions for a single scv face
     */
    PrimaryVariables evalNeumannSegment_(const SubControlVolumeFace &scvFace,
                                         const BoundaryTypes &bcTypes)
    {
        // temporary vector to store the neumann boundary fluxes
        PrimaryVariables flux(0);

        auto neumannFluxes = this->problem_().neumann(this->element_(), scvFace);

        // multiply neumann fluxes with the area and the extrusion factor
        const auto& scv = this->problem_().model().fvGeometries().subControlVolume(scvFace.insideScvIdx())
        neumannFluxes *= scvFace.area()*problem_().model().curVolVars(scv).extrusionFactor();

        // add fluxes to the residual
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            if (bcTypes.isNeumann(eqIdx))
            {
                int pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                flux[pvIdx] = neumannFluxes[eqIdx];
            }

        return flux;
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

    /*!
     * \brief Treat Dirichlet boundary conditions in a weak sense for a single
     *        intersection that only has Dirichlet boundary conditions
     */
    PrimaryVariables evalDirichletSegment_(const SubControlVolumeFace &scvFace,
                                           const BoundaryTypes &bcTypes)
    {
        // temporary vector to store the Dirichlet boundary fluxes
        PrimaryVariables flux(0);

        auto dirichletFlux = this->asImp_().computeFlux(scvFace);

        // multiply dirichlet fluxes with the area and the extrusion factor
        const auto& scv = this->problem_().model().fvGeometries().subControlVolume(scvFace.insideScvIdx())
        dirichletFlux *= scvFace.area()*problem_().model().curVolVars(scv).extrusionFactor();

        // add fluxes to the residual
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (bcTypes.isDirichlet(eqIdx))
            {
                int pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                flux[pvIdx] = dirichletFlux[eqIdx];
            }
        }

        return flux;
    }

    /*!
     * \brief Treat Dirichlet boundary conditions in a strong sense for a
     *        single intersection that has mixed D/N boundary conditions
     */
    void evalDirichletSegmentMixed_(PrimaryVariables &flux,
                                    const SubControlVolumeFace &scvFace)
    {
        // temporary vector to store the Dirichlet values
        PrimaryVariables dirichletValues = this->problem_().dirichlet(this->element_(), scvFace);

        // set Dirichlet conditions in a strong sense
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (bcTypes.isDirichlet(eqIdx))
            {
                int pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                this->residual_[0][eqIdx] = this->curPriVar_(0, pvIdx) - dirichletValues[pvIdx];
                flux[pvIdx] = 0;
            }
        }
    }

    /*!
     * \brief Add the flux terms to the local residual of the current element
     */
    void evalFluxes_()
    {
        // calculate the mass flux over the faces and subtract
        // it from the local rates
        int fIdx = -1;
        for (const auto& intersection : intersections(this->gridView_(), this->element_())) {
            if (!intersection.neighbor())
                continue;

            fIdx++;
            PrimaryVariables flux;

            Valgrind::SetUndefined(flux);
            this->asImp_().computeFlux(flux, fIdx);
            Valgrind::CheckDefined(flux);

            flux *= this->curVolVars_(0).extrusionFactor();

            this->residual_[0] += flux;
        }
    }
};

}

#endif   // DUMUX_CC_TPFA_LOCAL_RESIDUAL_HH