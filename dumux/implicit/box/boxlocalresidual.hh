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
#ifndef DUMUX_BOX_LOCAL_RESIDUAL_HH
#define DUMUX_BOX_LOCAL_RESIDUAL_HH

#include <dune/geometry/type.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/implicit/common/implicitlocalresidual.hh>

#include "boxproperties.hh"

namespace Dumux
{
/*!
 * \ingroup BoxModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the residual for models
 *        based on the fully implicit box scheme.
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class BoxLocalResidual : public ImplicitLocalResidual<TypeTag>
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

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef typename Dune::ReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<CoordScalar, dim> ReferenceElement;

    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    // copying the local residual class is not a good idea
    BoxLocalResidual(const BoxLocalResidual &);

public:
    BoxLocalResidual() : ParentType()
    { }

protected:
    /*!
     * \brief Set the values of the Dirichlet boundary control volumes
     *        of the current element.
     */
    void evalDirichlet_()
    {
        PrimaryVariables dirichletValues(0);
        for (int scvIdx = 0; scvIdx < this->fvGeometry_().numScv; ++scvIdx) {
            const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);
            
            if (bcTypes.hasDirichlet()) {
                // ask the problem for the dirichlet values
                const VertexPointer vPtr = this->element_().template subEntity<dim>(scvIdx);
                Valgrind::SetUndefined(dirichletValues);
                this->asImp_().problem_().dirichlet(dirichletValues, *vPtr);

                // set the dirichlet conditions
                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                    if (bcTypes.isDirichlet(eqIdx)) {
                        int pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                        assert(0 <= pvIdx && pvIdx < numEq);
                        Valgrind::CheckDefined(dirichletValues[pvIdx]);

                        this->residual_[scvIdx][eqIdx] =
                                this->curPriVar_(scvIdx, pvIdx) - dirichletValues[pvIdx];

                        this->storageTerm_[scvIdx][eqIdx] = 0.0;
                    }
                }
            }
        }
    }

    /*!
     * \brief Add all Neumann and outflow boundary conditions to the local
     *        residual.
     */
    void evalBoundaryFluxes_()
    {
        Dune::GeometryType geoType = this->element_().geometry().type();
        const ReferenceElement &refElement = ReferenceElements::general(geoType);

        IntersectionIterator isIt = this->gridView_().ibegin(this->element_());
        const IntersectionIterator &isEndIt = this->gridView_().iend(this->element_());
        for (; isIt != isEndIt; ++isIt)
        {
            // handle only faces on the boundary
            if (isIt->boundary()) {
                // Assemble the boundary for all vertices of the current
                // face
                int fIdx = isIt->indexInInside();
                int numFaceVerts = refElement.size(fIdx, 1, dim);
                for (int faceVertexIdx = 0;
                    faceVertexIdx < numFaceVerts;
                    ++faceVertexIdx)
                {
                    int scvIdx = refElement.subEntity(fIdx,
                                                        1,
                                                        faceVertexIdx,
                                                        dim);

                    int boundaryFaceIdx =
                        this->fvGeometry_().boundaryFaceIndex(fIdx, faceVertexIdx);

                    // add the residual of all vertices of the boundary
                    // segment
                    this->asImp_().evalNeumannSegment_(isIt,
                                                       scvIdx,
                                                       boundaryFaceIdx);
                    // evaluate the outflow conditions at the boundary face
                    this->asImp_().evalOutflowSegment_(isIt,
                                                       scvIdx,
                                                       boundaryFaceIdx);
                }
            }
        }
    }

    /*!
     * \brief Add Neumann boundary conditions for a single sub-control
     *        volume face to the local residual.
     */
    void evalNeumannSegment_(const IntersectionIterator &isIt,
                             const int scvIdx,
                             const int boundaryFaceIdx)
    {
        // temporary vector to store the neumann boundary fluxes
        PrimaryVariables neumannFlux(0.0);
        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);

        // deal with neumann boundaries
        if (bcTypes.hasNeumann()) {
            Valgrind::SetUndefined(neumannFlux);
            this->problem_().solDependentNeumann(neumannFlux,
                                          this->element_(),
                                          this->fvGeometry_(),
                                          *isIt,
                                          scvIdx,
                                          boundaryFaceIdx,
                                          this->curVolVars_());
            neumannFlux *=
                this->fvGeometry_().boundaryFace[boundaryFaceIdx].area
                * this->curVolVars_(scvIdx).extrusionFactor();
            Valgrind::CheckDefined(neumannFlux);

            // set the neumann conditions
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (!bcTypes.isNeumann(eqIdx))
                    continue;
                this->residual_[scvIdx][eqIdx] += neumannFlux[eqIdx];
            }
        }
    }

    /*!
    * \brief Add outflow boundary conditions for a single sub-control
    *        volume face to the local residual.
    *
    * \param isIt   The intersection iterator of current element
    * \param scvIdx The index of the considered face of the sub-control volume
    * \param boundaryFaceIdx The index of the considered boundary face of the sub control volume
    */
    void evalOutflowSegment_(const IntersectionIterator &isIt,
                            const int scvIdx,
                            const int boundaryFaceIdx)
    {
        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);
        // deal with outflow boundaries
        if (bcTypes.hasOutflow())
        {
            //calculate outflow fluxes
            PrimaryVariables values(0.0);
            this->asImp_().computeFlux(values, boundaryFaceIdx, true);
            Valgrind::CheckDefined(values);
            
            for (int equationIdx = 0; equationIdx < numEq; ++equationIdx)
            {
                if (!bcTypes.isOutflow(equationIdx) )
                    continue;
                // deduce outflow
                this->residual_[scvIdx][equationIdx] += values[equationIdx];
            }
        }
    }

    /*!
     * \brief Add the flux terms to the local residual of all
     *        sub-control volumes of the current element.
     */
    void evalFluxes_()
    {
        // calculate the mass flux over the faces and subtract
        // it from the local rates
        for (int scvfIdx = 0; scvfIdx < this->fvGeometry_().numScvf; scvfIdx++)
        {
            int i = this->fvGeometry_().subContVolFace[scvfIdx].i;
            int j = this->fvGeometry_().subContVolFace[scvfIdx].j;

            PrimaryVariables flux;

            Valgrind::SetUndefined(flux);
            this->asImp_().computeFlux(flux, scvfIdx);
            Valgrind::CheckDefined(flux);

            Scalar extrusionFactor =
                (this->curVolVars_(i).extrusionFactor()
                 + this->curVolVars_(j).extrusionFactor())
                / 2;
            flux *= extrusionFactor;

            // The balance equation for a finite volume is:
            //
            // dStorage/dt = Flux + Source
            //
            // where the 'Flux' and the 'Source' terms represent the
            // mass per second which _ENTER_ the finite
            // volume. Re-arranging this, we get
            //
            // dStorage/dt - Source - Flux = 0
            //
            // Since the flux calculated by computeFlux() goes _OUT_
            // of sub-control volume i and _INTO_ sub-control volume
            // j, we need to add the flux to finite volume i and
            // subtract it from finite volume j
            this->residual_[i] += flux;
            this->residual_[j] -= flux;
        }
    }
    
    const VertexMapper &vertexMapper_() const
    { return this->problem_().vertexMapper(); };
};

}

#endif
