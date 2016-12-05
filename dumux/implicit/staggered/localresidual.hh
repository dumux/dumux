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
#ifndef DUMUX_STAGGERED_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>

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
class StaggeredLocalResidual
{
    using ParentType = ImplicitLocalResidual<TypeTag>;
    friend class ImplicitLocalResidual<TypeTag>;
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
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);
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

public:
    // copying the local residual class is not a good idea
    StaggeredLocalResidual(const StaggeredLocalResidual &) = delete;

    StaggeredLocalResidual() = default;


     /*!
     * \brief Initialize the local residual.
     *
     * This assumes that all objects of the simulation have been fully
     * allocated but not necessarily initialized completely.
     *
     * \param problem The representation of the physical problem to be
     *             solved.
     */
    void init(Problem &problem)
    { problemPtr_ = &problem; }


     /*!
     * \name User interface
     * \note The following methods are usually expensive to evaluate
     *       They are useful for outputting residual information.
     */
    // \{

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     */
    void eval(const Element &element)
    {
        // make sure FVElementGeometry and volume variables are bound to the element
        auto fvGeometry = localView(this->problem().model().globalFvGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(problem().model().curGlobalVolVars());
        curElemVolVars.bind(element, fvGeometry, problem().model().curSol());

        auto prevElemVolVars = localView(problem().model().prevGlobalVolVars());
        prevElemVolVars.bindElement(element, fvGeometry, problem().model().prevSol());

        auto elemFluxVarsCache = localView(problem().model().globalFluxVarsCache());
        elemFluxVarsCache.bindElement(element, fvGeometry, curElemVolVars);

        ElementBoundaryTypes bcTypes;
        bcTypes.update(problem(), element, fvGeometry);

        auto& curGlobalFaceVars = problem().model().curGlobalFaceVars();
        auto& prevGlobalFaceVars = problem().model().prevGlobalFaceVars();


        asImp_().eval(element, fvGeometry,
                      prevElemVolVars, curElemVolVars,
                      prevGlobalFaceVars, curGlobalFaceVars,
                      bcTypes, elemFluxVarsCache);
    }

     /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        equations from zero.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevVolVars The volume averaged variables for all
     *                   sub-control volumes of the element at the previous
     *                   time level
     * \param curVolVars The volume averaged variables for all
     *                   sub-control volumes of the element at the current
     *                   time level
     * \param bcTypes The types of the boundary conditions for all
     *                vertices of the element
     */
    void eval(const Element &element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& prevElemVolVars,
              const ElementVolumeVariables& curElemVolVars,
              const GlobalFaceVars& prevFaceVars,
              const GlobalFaceVars& curFaceVars,
              const ElementBoundaryTypes &bcTypes,
              const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // resize the vectors for all terms
        const auto numScv = fvGeometry.numScv();
        const auto numScvf = fvGeometry.numScvf();
        ccResidual_.resize(numScv, 0.0);
        ccStorageTerm_.resize(numScv, 0.0);

        faceResiduals_.resize(numScvf);
        faceStorageTerms_.resize(numScvf);
        for(int i = 0; i < numScvf; ++i)
        {
            faceResiduals_[i].resize(1, 0.0);
            faceStorageTerms_[i].resize(1, 0.0);
        }

        asImp_().evalVolumeTerms_(element, fvGeometry, prevElemVolVars, curElemVolVars, prevFaceVars, curFaceVars, bcTypes);
        asImp_().evalFluxes_(element, fvGeometry, curElemVolVars, curFaceVars, bcTypes, elemFluxVarsCache);
        asImp_().evalBoundary_(element, fvGeometry, curElemVolVars, bcTypes, elemFluxVarsCache);
    }

     /*!
     * \brief Return the problem we are solving. Only call this after init()!
     */
    const Problem& problem() const
    { return *problemPtr_; }

    /*!
     * \brief Return the problem we are solving. Only call this after init()!
     */
    Problem& problem()
    { return *problemPtr_; }


protected:

    void evalFluxes_(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const GlobalFaceVars& faceVars,
                     const ElementBoundaryTypes& bcTypes,
                     const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        evalFluxesForCellCenter_(element, fvGeometry, elemVolVars, faceVars, bcTypes, elemFluxVarsCache);
        evalFluxesForFaces_(element, fvGeometry, elemVolVars, faceVars, bcTypes, elemFluxVarsCache);
    }

    void evalFluxesForCellCenter_(const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& elemVolVars,
                                  const GlobalFaceVars& faceVars,
                                  const ElementBoundaryTypes& bcTypes,
                                  const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        for (auto&& scvf : scvfs(fvGeometry))
        {
            ccResidual_[0] += asImp_().computeFluxForCellCenter_(element, fvGeometry, elemVolVars, faceVars, scvf, bcTypes, elemFluxVarsCache[scvf]);
        }
    }

    void evalFluxesForFaces_(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const GlobalFaceVars& globalFaceVars,
                             const ElementBoundaryTypes& bcTypes,
                             const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if(!scvf.boundary())
                faceResiduals_[scvf.localFaceIdx()][0] += asImp_().computeFluxForFace(scvf, elemVolVars, globalFaceVars);
        }
    }



    CellCenterPrimaryVariables computeFluxForCellCenter_(const Element &element,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVolumeVariables& elemVolVars,
                                               const GlobalFaceVars& faceVars,
                                               const SubControlVolumeFace &scvf,
                                               const ElementBoundaryTypes& bcTypes,
                                               const FluxVariablesCache& fluxVarsCache)
    {
//         if (!scvf.boundary() /*TODO: || GET_PROP_VALUE(TypeTag, BoundaryReconstruction)*/)
//             return /*this->asImp_().*/asImp_().computeFluxForCellCenter(element, fvGeometry, elemVolVars, faceVars, scvf, fluxVarsCache);
//         else
//         {
//             return CellCenterPrimaryVariables(0.0);
//         }

         return asImp_().computeFluxForCellCenter(element, fvGeometry, elemVolVars, faceVars, scvf, fluxVarsCache);
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
//                 auto test = this->problem().faceDirichletAtPos(scvf.center(), scvf.directionIndex());
//              std::cout << " face: " << scvf.center() << ",  normal: "  <<  scvf.unitOuterNormal()    << ",   value: " << test << std::endl;

//                 auto bcTypes = this->problem().boundaryTypes(element, scvf);

//                 this->residual_[0] += evalBoundaryFluxes_(element, fvGeometry, elemVolVars, scvf, bcTypes, elemFluxVarsCache[scvf]);
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

        auto neumannFluxes = this->problem().neumann(element, fvGeometry, elemVolVars, scvf);

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
//         PrimaryVariables dirichletValues = this->problem().dirichlet(element, scvf);
//
//         // get the primary variables
//         const auto& priVars = elemVolVars[scvf.insideScvIdx()].priVars();
//
//         // set Dirichlet conditions in a strong sense
//         for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
//         {
//             if (bcTypes.isDirichlet(eqIdx))
//             {
//                 auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
//                 this->residual_[0][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];
//             }
//         }
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
     * \brief Add the change the source term for stationary problems
     *        to the local residual of all sub-control volumes of the
     *        current element.
     */
    template<class P = Problem>
    typename std::enable_if<Dumux::Capabilities::isStationary<P>::value, void>::type
    evalVolumeTerms_(const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const GlobalFaceVars& prevFaceVars,
                     const GlobalFaceVars& curFaceVars,
                     const ElementBoundaryTypes &bcTypes)
    {
        // evaluate the volume terms (storage + source terms)
        for (auto&& scv : scvs(fvGeometry))
        {
            auto curExtrusionFactor = curElemVolVars[scv].extrusionFactor();

            // subtract the source term from the local rate
            CellCenterPrimaryVariables source = asImp_().computeSourceForCellCenter(element, fvGeometry, curElemVolVars, curFaceVars, scv);
            source *= scv.volume()*curExtrusionFactor;

            ccResidual_[0] -= source;

            // now, treat the dofs on the facets:
            for(auto&& scvf : scvfs(fvGeometry))
            {
                // the source term:
                auto faceSource = asImp_().computeSourceForFace(scvf, curElemVolVars, curFaceVars);
                faceSource *= 0.5*scv.volume()*curExtrusionFactor;
                faceResiduals_[scvf.localFaceIdx()][0] -= faceSource;
            }
        }
    }

    /*!
     * \brief Add the change in the storage terms and the source term
     *        to the local residual of all sub-control volumes of the
     *        current element.
     */
    template<class P = Problem>
    typename std::enable_if<!Dumux::Capabilities::isStationary<P>::value, void>::type
    evalVolumeTerms_(const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& prevElemVolVars,
                     const ElementVolumeVariables& curElemVolVars,
                     const GlobalFaceVars& prevFaceVars,
                     const GlobalFaceVars& curFaceVars,
                     const ElementBoundaryTypes &bcTypes)
    {
        // evaluate the volume terms (storage + source terms)
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& curVolVars = curElemVolVars[scv];
            const auto& prevVolVars = prevElemVolVars[scv];

            // mass balance within the element. this is the
            // \f$\frac{m}{\partial t}\f$ term if using implicit
            // euler as time discretization.
            //
            // We might need a more explicit way for
            // doing the time discretization...
            auto prevCCStorage = asImp_().computeStorageForCellCenter(scv, prevVolVars);
            auto curCCStorage = asImp_().computeStorageForCellCenter(scv, curVolVars);

            prevCCStorage *= prevVolVars.extrusionFactor();
            curCCStorage *= curVolVars.extrusionFactor();

            ccStorageTerm_[0] = std::move(curCCStorage);
            ccStorageTerm_[0] -= std::move(prevCCStorage);
            ccStorageTerm_[0] *= scv.volume();
            ccStorageTerm_[0] /= problem().timeManager().timeStepSize();

            // add the storage term to the residual
            ccResidual_[0] += ccStorageTerm_[0];

            // subtract the source term from the local rate
            CellCenterPrimaryVariables source = asImp_().computeSourceForCellCenter(element, fvGeometry, curElemVolVars, curFaceVars, scv);
            source *= scv.volume()*curVolVars.extrusionFactor();

            ccResidual_[0] -= source;


            // now, treat the dofs on the facets:
            for(auto&& scvf : scvfs(fvGeometry))
            {
                // the storage term:
                auto prevFaceStorage = asImp_().computeStorageForFace(scvf, prevVolVars, prevFaceVars);
                auto curFaceStorage = asImp_().computeStorageForFace(scvf, curVolVars, prevFaceVars);
                faceStorageTerms_[scvf.localFaceIdx()][0] = std::move(curFaceStorage);
                faceStorageTerms_[scvf.localFaceIdx()][0] -= std::move(prevFaceStorage);
                faceStorageTerms_[scvf.localFaceIdx()][0] *= (scv.volume()/2.0);
                faceStorageTerms_[scvf.localFaceIdx()][0] /= problem().timeManager().timeStepSize();

                faceResiduals_[scvf.localFaceIdx()][0] += faceStorageTerms_[scvf.localFaceIdx()][0];

                // the source term:
                auto faceSource = asImp_().computeSourceForFace(scvf, curElemVolVars, curFaceVars);
                faceSource *= 0.5*scv.volume()*curVolVars.extrusionFactor();
                faceResiduals_[scvf.localFaceIdx()][0] -= faceSource;
            }
        }
    }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    CellCenterSolutionVector ccResidual_;
    CellCenterSolutionVector ccStorageTerm_;
    std::vector<FaceSolutionVector> faceResiduals_;
    std::vector<FaceSolutionVector> faceStorageTerms_;

private:
    Problem* problemPtr_;

};

}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH
