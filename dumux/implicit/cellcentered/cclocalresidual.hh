// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \brief Calculates the residual of models based on the box scheme element-wise.
 */
#ifndef DUMUX_CC_LOCAL_RESIDUAL_HH
#define DUMUX_CC_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>
#include <dune/grid/common/geometry.hh>

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
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    // copying the local residual class is not a good idea
    CCLocalResidual(const CCLocalResidual &);

public:
    CCLocalResidual() : ParentType()
    { }

protected:

    /*!
     * \brief Set the values of the Dirichlet boundary control volumes
     *        of the current element.
     */
    void evalDirichlet_()
    {
        PrimaryVariables dirichletValues(0);

        BoundaryTypes bcTypes;
        IntersectionIterator isIt = this->problem_().gridView().ibegin(this->element_());
        IntersectionIterator isEndIt = this->problem_().gridView().iend(this->element_());
        for (; isIt != isEndIt; ++isIt) {
            if (!isIt->boundary())
                continue;
            
            this->problem_().boundaryTypes(bcTypes, *isIt);
            if (bcTypes.hasDirichlet())
                break;
        }
        
        Valgrind::SetUndefined(dirichletValues);
        // HACK: ask for Dirichlet value at element center
        this->asImp_().problem_().dirichletAtPos(dirichletValues, this->element_().geometry().center());

        // set the dirichlet conditions
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!bcTypes.isDirichlet(eqIdx))
                continue;
            int pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
            assert(0 <= pvIdx && pvIdx < numEq);
            Valgrind::CheckDefined(dirichletValues[pvIdx]);

            this->residual_[0][eqIdx] =
                    (this->curPriVar_(0, pvIdx) - dirichletValues[pvIdx]);

            this->storageTerm_[0][eqIdx] = 0.0;
        }
    }

    /*!
     * \brief Add all Neumann and outflow boundary conditions to the local
     *        residual.
     */
    void evalBoundaryFluxes_()
    {
        IntersectionIterator isIt = this->gridView_().ibegin(this->element_());
        const IntersectionIterator &endIt = this->gridView_().iend(this->element_());
        for (; isIt != endIt; ++isIt)
        {
            // handle only faces on the boundary
            if (!isIt->boundary())
                continue;

            // add the residual of all vertices of the boundary
            // segment
            this->asImp_().evalNeumannSegment_(isIt);

            // evaluate the outflow conditions at the boundary face
            this->asImp_().evalOutflowSegment_(isIt);
        }
    }

    /*!
     * \brief Add Neumann boundary conditions for a single intersection
     */
    void evalNeumannSegment_(const IntersectionIterator &isIt)
    {
        // temporary vector to store the neumann boundary fluxes
        PrimaryVariables values(0.0);

        BoundaryTypes bcTypes;
        this->problem_().boundaryTypes(bcTypes, *isIt);

        // deal with neumann boundaries
        if (bcTypes.hasNeumann()) {
            Valgrind::SetUndefined(values);
            this->problem_().boxSDNeumann(values,
                                          this->element_(),
                                          this->fvGeometry_(),
                                          *isIt,
                                          /*scvIdx=*/0,
                                          /*boundaryFaceIdx=*/isIt->indexInInside(),
                                          this->curVolVars_());
            values *= isIt->geometry().volume()
                * this->curVolVars_(0).extrusionFactor();
            Valgrind::CheckDefined(values);

            // set the neumann conditions
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (!bcTypes.isNeumann(eqIdx))
                    continue;
                this->residual_[0][eqIdx] += values[eqIdx];
            }
        }
    }

    /*!
     * \brief Add outflow boundary conditions for a single intersection
     */
    void evalOutflowSegment_(const IntersectionIterator &isIt)
    {
        BoundaryTypes bcTypes;
        this->problem_().boundaryTypes(bcTypes, *isIt);
        
        // deal with outflow boundaries
        if (bcTypes.hasOutflow())
        {
            PrimaryVariables values(0.0);
            this->asImp_().computeFlux(values, 
                                       /*boundaryFaceIdx=*/isIt->indexInInside(), 
                                       true);
            Valgrind::CheckDefined(values);
            
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                if (!bcTypes.isOutflow(eqIdx) )
                    continue;
                this->residual_[0][eqIdx] += values[eqIdx];
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
        int faceIdx = -1;
        IntersectionIterator isIt = this->gridView_().ibegin(this->element_());
        IntersectionIterator isEndIt = this->gridView_().iend(this->element_());
        for (; isIt != isEndIt; ++isIt) {
            if (!isIt->neighbor())
                continue;

            faceIdx++;
	    PrimaryVariables flux;

            Valgrind::SetUndefined(flux);
            this->asImp_().computeFlux(flux, faceIdx);
            Valgrind::CheckDefined(flux);

            flux *= this->curVolVars_(0).extrusionFactor();

            this->residual_[0] += flux;
        }
    }
};

}

#endif
