/*****************************************************************************
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
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
 * \brief Caculates the residual of models based on the box scheme element-wise.
 */
#ifndef DUMUX_BOX_LOCAL_RESIDUAL_HH
#define DUMUX_BOX_LOCAL_RESIDUAL_HH

#include <dune/istl/bvector.hh>
#include <dune/grid/common/geometry.hh>

#include <dumux/common/valgrind.hh>

#include "boxproperties.hh"

namespace Dumux
{
/*!
 * \ingroup BoxModel
 * \brief Element-wise caculation of the residual matrix for models
 *        based on the box scheme .
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class BoxLocalResidual
{
private:
    typedef BoxLocalResidual<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalResidual)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum { dim = GridView::dimension };
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename GridView::ctype CoordScalar;

    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;

    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Element::Geometry Geometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVariables)) ElementVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    enum { numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)) };

    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> LocalBlockVector;
    
    // copying the local residual class is not a good idea
    BoxLocalResidual(const BoxLocalResidual &)
    {};

public:
    BoxLocalResidual()
    { }

    ~BoxLocalResidual()
    { }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        conservation equations from zero and store the results
     *        internally.
     *
     * The results can be requested afterwards using the residual()
     * and storageTerm() methods.
     */
    void eval(const ElementVariables &elemVars)
    {
        int numScv = elemVars.numScv();
        internalResidual_.resize(numScv);
        internalStorageTerm_.resize(numScv);
        eval(internalResidual_, internalStorageTerm_, elemVars);
    };
    
    /*!
     * \brief Return the result of the eval() call using internal
     *        storage.
     */
    const LocalBlockVector &residual() const
    { return internalResidual_; }
    const VectorBlock &residual(int scvIdx) const
    { return internalResidual_[scvIdx]; }

    /*!
     * \brief Return the storage term calculated using the last call
     *        to eval() using internal storage.
     */
    const LocalBlockVector &storageTerm() const
    { return internalStorageTerm_; }
    const VectorBlock &storageTerm(int scvIdx) const
    { return internalStorageTerm_[scvIdx]; }
    
    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        conservation equations from zero.
     *
     * \param residual Stores the residual vector
     * \param storageTerm Stores the derivative of the storage term to time
     * \param elemVars All the element's secondary variables required to calculate the local residual
     */
    void eval(LocalBlockVector &residual,
              LocalBlockVector &storageTerm,
              const ElementVariables &elemVars) const
    {
        residual = 0.0;
        storageTerm = 0.0;

        // evaluate the flux terms
        asImp_().evalFluxes(residual, elemVars);

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < elemVars.fvElemGeom().numVertices; i++)
            Valgrind::CheckDefined(residual[i]);
#endif // HAVE_VALGRIND

        // evaluate the storage and the source terms
        asImp_().evalVolumeTerms_(residual, storageTerm, elemVars);

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < elemVars.fvElemGeom().numVertices; i++) 
            Valgrind::CheckDefined(residual[i]);
#endif // !defined NDEBUG && HAVE_VALGRIND

        // evaluate the boundary conditions
        asImp_().evalBoundary_(residual, storageTerm, elemVars);

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < elemVars.fvElemGeom().numVertices; i++)
            Valgrind::CheckDefined(residual[i]);
#endif // HAVE_VALGRIND
    }

    /*!
     * \brief Calculate the amount of all conservation quantities
     *        stored in all element's sub-control volumes for a given
     *        history index.
     *
     * This is used to figure out how much of each conservation
     * quantity is inside the element.
     */
    void evalStorage(const ElementVariables &elemVars, 
                     int historyIdx)
    {
        int numScv = elemVars.numScv();
        internalStorageTerm_.resize(numScv);
        evalStorage(internalStorageTerm_,
                    elemVars,
                    historyIdx);
    }

    /*!
     * \brief Calculate the amount of all conservation quantities
     *        stored in all element's sub-control volumes for a given
     *        history index.
     *
     * This is used to figure out how much of each conservation
     * quantity is inside the element.
     */
    void evalStorage(LocalBlockVector &storage,
                     const ElementVariables &elemVars, 
                     int historyIdx) const
    {
        // calculate the amount of conservation each quantity inside
        // all sub control volumes
        for (int scvIdx=0; scvIdx < elemVars.numScv(); scvIdx++)
        {
            storage[scvIdx] = 0.0;
            asImp_().computeStorage(storage[scvIdx], 
                                    elemVars,
                                    scvIdx,
                                    historyIdx);
            storage[scvIdx] *= 
                elemVars.fvElemGeom().subContVol[scvIdx].volume
                * elemVars.volVars(scvIdx).extrusionFactor();
        }
    }

    /*!
     * \brief Add the flux term to a local residual.
     */
    void evalFluxes(LocalBlockVector &residual,
                    const ElementVariables &elemVars) const
    {
        PrimaryVariables flux;
        
        // calculate the mass flux over the sub-control volume faces
        for (int scvfIdx = 0; 
             scvfIdx < elemVars.fvElemGeom().numEdges;
             scvfIdx++)
        {
            int i = elemVars.fvElemGeom().subContVolFace[scvfIdx].i;
            int j = elemVars.fvElemGeom().subContVolFace[scvfIdx].j;

            Valgrind::SetUndefined(flux);
            asImp_().computeFlux(flux, /*context=*/elemVars, scvfIdx);
            flux *= elemVars.fluxVars(scvfIdx).extrusionFactor();
            Valgrind::CheckDefined(flux);

            // The balance equation for a finite volume is given by
            //
            // dStorage/dt = Flux + Source
            //
            // where the 'Flux' and the 'Source' terms represent the
            // mass per second which _ENTER_ the finite
            // volume. Re-arranging this, we get
            //
            // dStorage/dt - Source - Flux = 0
            //
            // Since the mass flux as calculated by computeFlux() goes
            // _OUT_ of sub-control volume i and _INTO_ sub-control
            // volume j, we need to add the flux to finite volume i
            // and subtract it from finite volume j
            residual[i] += flux;
            residual[j] -= flux;
        }
    }

protected:
    /*!
     * \brief Evaluate the boundary conditions of an element.
     */
    void evalBoundary_(LocalBlockVector &residual,
                       LocalBlockVector &storageTerm,
                       const ElementVariables &elemVars) const
    {
        if (!elemVars.onBoundary()) {
            return;
        }
        
        if (elemVars.hasNeumann())
            asImp_().evalNeumann_(residual, elemVars);
        
        if (elemVars.hasDirichlet())
            asImp_().evalDirichlet_(residual, storageTerm, elemVars);
    }

    /*!
     * \brief Set the values of the Dirichlet boundary control volumes
     *        of the current element.
     */
    void evalDirichlet_(LocalBlockVector &residual,
                        LocalBlockVector &storageTerm,
                        const ElementVariables &elemVars) const
    {
        PrimaryVariables tmp(0);
        for (int scvIdx = 0; scvIdx < elemVars.numScv(); ++scvIdx) {
            const BoundaryTypes &bcTypes = elemVars.boundaryTypes(scvIdx);
            if (!bcTypes.hasDirichlet())
                continue;
            
            // ask the problem for the dirichlet values
            asImp_().computeDirichlet_(tmp, elemVars, scvIdx);

            const PrimaryVariables &priVars = 
                elemVars.volVars(scvIdx).primaryVars();

            // set the dirichlet conditions
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (!bcTypes.isDirichlet(eqIdx))
                    continue;
                int pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                
                assert(0 <= pvIdx && pvIdx < numEq);
                Valgrind::CheckDefined(pvIdx);
                Valgrind::CheckDefined(priVars);
                Valgrind::CheckDefined(tmp[pvIdx]);

                residual[scvIdx][eqIdx] = priVars[pvIdx] - tmp[pvIdx];
                storageTerm[scvIdx][eqIdx] = 0.0;
            };
        };
    }

    /*!
     * \brief Add all Neumann boundary conditions to the local
     *        residual.
     */
    void evalNeumann_(LocalBlockVector &residual,
                      const ElementVariables &elemVars) const
    {
        const Element &elem = elemVars.element();
        Dune::GeometryType geoType = elem.geometry().type();
        const ReferenceElement &refElem = ReferenceElements::general(geoType);

        const GridView &gridView = elemVars.gridView();
        IntersectionIterator isIt = gridView.ibegin(elem);
        const IntersectionIterator &endIt = gridView.iend(elem);
        for (; isIt != endIt; ++isIt)
        {
            // handle only faces on the boundary
            if (!isIt->boundary())
                continue;

            // Assemble the boundary for all vertices of the current
            // face
            int faceIdx = isIt->indexInInside();
            int numFaceVerts = refElem.size(faceIdx, 1, dim);
            for (int faceVertIdx = 0;
                 faceVertIdx < numFaceVerts;
                 ++faceVertIdx)
            {
                int scvIdx = refElem.subEntity(/*entityIdx=*/faceIdx,
                                               /*entityCodim=*/1,
                                               /*subEntityIdx=*/faceVertIdx,
                                               /*subEntityCodim=*/dim);
                
                int boundaryFaceIdx =
                    elemVars.fvElemGeom().boundaryFaceIndex(faceIdx, faceVertIdx);

                // add the residual of all vertices of the boundary
                // segment
                evalNeumannSegment_(residual,
                                    elemVars,
                                    *isIt,
                                    scvIdx,
                                    boundaryFaceIdx);
            }
        }
    }

    void computeDirichlet_(PrimaryVariables &values, 
                           const ElementVariables &elemVars,
                           int scvIdx) const
    { elemVars.problem().dirichlet(values, elemVars, scvIdx); }


    /*!
     * \brief Add Neumann boundary conditions for a single sub-control
     *        volume face to the local residual.
     */
    void evalNeumannSegment_(LocalBlockVector &residual,
                             const ElementVariables &elemVars,
                             const Intersection &is,
                             int scvIdx,
                             int boundaryFaceIdx) const
    {
        // temporary vector to store the neumann boundary fluxes
        const BoundaryTypes &bcTypes = elemVars.boundaryTypes(scvIdx);
        PrimaryVariables values;
        
        // deal with neumann boundaries
        if (bcTypes.hasNeumann()) {
            elemVars.problem().neumann(values,
                                       elemVars,
                                       is,
                                       scvIdx,
                                       boundaryFaceIdx);
            Valgrind::CheckDefined(values);

            values *= 
                elemVars.fvElemGeom().boundaryFace[boundaryFaceIdx].area
                * elemVars.volVars(scvIdx, /*historyIdx*/0).extrusionFactor();
            Valgrind::CheckDefined(values);
            residual[scvIdx] += values;
        }
    }

    /*!
     * \brief Add the change in the storage terms and the source term
     *        to the local residual of all sub-control volumes of the
     *        current element.
     */
    void evalVolumeTerms_(LocalBlockVector &residual, 
                          LocalBlockVector &storageTerm,
                          const ElementVariables &elemVars) const
    {
        PrimaryVariables tmp, tmp2;

        // evaluate the volume terms (storage + source terms)
        for (int scvIdx=0; scvIdx < elemVars.numScv(); scvIdx++)
        {
            Scalar extrusionFactor =
                elemVars.volVars(scvIdx, /*historyIdx=*/0).extrusionFactor();          

            // mass balance within the element. this is the
            // $\frac{m}{\partial t}$ term if using implicit
            // euler as time discretization.
            //
            // TODO (?): we might need a more explicit way for
            // doing the time discretization...
            asImp_().computeStorage(tmp2, 
                                    elemVars,
                                    scvIdx, 
                                    /*historyIdx=*/1);
            asImp_().computeStorage(tmp, 
                                    elemVars,
                                    scvIdx,
                                    /*historyIdx=*/0);

            tmp -= tmp2;
            tmp *=
                elemVars.fvElemGeom().subContVol[scvIdx].volume
                / elemVars.problem().timeManager().timeStepSize()
                * extrusionFactor;

            storageTerm[scvIdx] += tmp;
            residual[scvIdx] += tmp;
            
            // subtract the source term from the residual
            asImp_().computeSource(tmp2, elemVars, scvIdx);
            tmp2 *= elemVars.fvElemGeom().subContVol[scvIdx].volume * extrusionFactor;
            residual[scvIdx] -= tmp2;

            // make sure that only defined quantities were used
            // to calculate the residual.
            Valgrind::CheckDefined(storageTerm[scvIdx]);
            Valgrind::CheckDefined(residual[scvIdx]);
        }
    }

private:
    Implementation &asImp_()
    {
      assert(static_cast<Implementation*>(this) != 0);
      return *static_cast<Implementation*>(this);
    }

    const Implementation &asImp_() const
    {
      assert(static_cast<const Implementation*>(this) != 0);
      return *static_cast<const Implementation*>(this);
    }

    LocalBlockVector internalResidual_;
    LocalBlockVector internalStorageTerm_;
};

}

#endif
