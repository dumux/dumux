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
#include <dumux/implicit/common/implicitlocalresidual.hh>

#include "ccproperties.hh"

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
    typedef ImplicitLocalResidual<TypeTag> ParentType;
    friend class ImplicitLocalResidual<TypeTag>;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq)
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    // copying the local residual class is not a good idea
    CCLocalResidual(const CCLocalResidual &);

public:
    CCLocalResidual() : ParentType()
    { }

protected:

    /*!
     * \brief Add all fluxes resulting from Neumann, outflow and pure Dirichlet
     *        boundary conditions to the local residual.
     */
    void evalBoundaryFluxes_()
    {
        IntersectionIterator isIt = this->gridView_().ibegin(this->element_());
        const IntersectionIterator &isEndIt = this->gridView_().iend(this->element_());
        for (; isIt != isEndIt; ++isIt)
        {
            // handle only faces on the boundary
            if (!isIt->boundary())
                continue;

            BoundaryTypes bcTypes;
            this->problem_().boundaryTypes(bcTypes, *isIt);

            // evaluate the Neumann conditions at the boundary face
            if (bcTypes.hasNeumann()) 
                this->asImp_().evalNeumannSegment_(isIt, bcTypes);

            // evaluate the outflow conditions at the boundary face
            if (bcTypes.hasOutflow()) 
                this->asImp_().evalOutflowSegment_(isIt, bcTypes);

            // evaluate the pure Dirichlet conditions at the boundary face
            if (bcTypes.hasDirichlet() && !bcTypes.hasNeumann()) 
                this->asImp_().evalDirichletSegmentPure_(isIt, bcTypes);
        }
    }

    /*!
     * \brief Evaluate Dirichlet conditions on faces that have mixed
     *        Dirichlet/Neumann boundary conditions.
     */
    void evalDirichlet_()
    {
        IntersectionIterator isIt = this->gridView_().ibegin(this->element_());
        const IntersectionIterator &isEndIt = this->gridView_().iend(this->element_());
        for (; isIt != isEndIt; ++isIt)
        {
            // handle only faces on the boundary
            if (!isIt->boundary())
                continue;

            BoundaryTypes bcTypes;
            this->problem_().boundaryTypes(bcTypes, *isIt);

            if (bcTypes.hasDirichlet() && bcTypes.hasNeumann()) 
                this->asImp_().evalDirichletSegmentMixed_(isIt, bcTypes);
        }
    }

    /*!
     * \brief Add Neumann boundary conditions for a single intersection
     */
    void evalNeumannSegment_(const IntersectionIterator &isIt, 
                             const BoundaryTypes &bcTypes)
    {
        // temporary vector to store the neumann boundary fluxes
        PrimaryVariables values;
        Valgrind::SetUndefined(values);

        unsigned bfIdx = isIt->indexInInside();
        this->problem_().solDependentNeumann(values,
                                      this->element_(),
                                      this->fvGeometry_(),
                                      *isIt,
                                      /*scvIdx=*/0,
                                      bfIdx,
                                      this->curVolVars_());
        values *= isIt->geometry().volume()
                 * this->curVolVars_(0).extrusionFactor();

        // add fluxes to the residual
        Valgrind::CheckDefined(values);
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (bcTypes.isNeumann(eqIdx))
                this->residual_[0][eqIdx] += values[eqIdx];
        }
    }

    /*!
     * \brief Add outflow boundary conditions for a single intersection
     */
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
    }

    /*!
     * \brief Treat Dirichlet boundary conditions in a weak sense for a single
     *        intersection that only has Dirichlet boundary conditions
     */
    void evalDirichletSegmentPure_(const IntersectionIterator &isIt, 
                                   const BoundaryTypes &bcTypes)
    {
        // temporary vector to store the Dirichlet boundary fluxes
        PrimaryVariables values;
        Valgrind::SetUndefined(values);

        unsigned bfIdx = isIt->indexInInside();

        this->asImp_().computeFlux(values, bfIdx, true);
        values *= this->curVolVars_(0).extrusionFactor();

        // add fluxes to the residual
        Valgrind::CheckDefined(values);
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (bcTypes.isDirichlet(eqIdx))
                this->residual_[0][eqIdx] += values[eqIdx];
        }
    }

    /*!
     * \brief Treat Dirichlet boundary conditions in a strong sense for a
     *        single intersection that has mixed D/N boundary conditions
     */
    void evalDirichletSegmentMixed_(const IntersectionIterator &isIt, 
                                    const BoundaryTypes &bcTypes)
    {
        // temporary vector to store the Dirichlet boundary fluxes
        PrimaryVariables values;
        Valgrind::SetUndefined(values);

        unsigned bfIdx = isIt->indexInInside();

        this->problem_().dirichlet(values, *isIt);
        Valgrind::CheckDefined(values);

        // set Dirichlet conditions in a strong sense
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (bcTypes.isDirichlet(eqIdx))
            {
                int pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                this->residual_[0][eqIdx] 
                  = this->curPriVar_(0, pvIdx) - values[pvIdx];
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
        IntersectionIterator isIt = this->gridView_().ibegin(this->element_());
        IntersectionIterator isEndIt = this->gridView_().iend(this->element_());
        for (; isIt != isEndIt; ++isIt) {
            if (!isIt->neighbor())
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

#endif
