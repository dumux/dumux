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
        // resize and reset all terms
        const auto numScvf = fvGeometry.numScvf();

        ccResidual_ = 0.0;
        ccStorageTerm_ = 0.0;

        faceResiduals_.resize(numScvf);
        faceStorageTerms_.resize(numScvf);
        faceResiduals_ = 0.0;
        faceStorageTerms_ = 0.0;

        asImp_().evalVolumeTerms_(element, fvGeometry, prevElemVolVars, curElemVolVars, prevFaceVars, curFaceVars, bcTypes);
        asImp_().evalFluxes_(element, fvGeometry, curElemVolVars, curFaceVars, bcTypes, elemFluxVarsCache);
        asImp_().evalBoundary_(element, fvGeometry, curElemVolVars, curFaceVars, bcTypes, elemFluxVarsCache);
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

    const auto& ccResidual() const
    { return ccResidual_; }

    const auto& faceResiduals() const
    { return faceResiduals_; }

    const auto& faceResidual(const int faceIdx) const
    { return faceResiduals_[faceIdx]; }


protected:

     /*!
     * \brief Evaluate the flux terms for both cell center and face dofs
     */
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

     /*!
     * \brief Evaluate the flux terms for cell center dofs
     */
    void evalFluxesForCellCenter_(const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& elemVolVars,
                                  const GlobalFaceVars& faceVars,
                                  const ElementBoundaryTypes& bcTypes,
                                  const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if(!scvf.boundary())
                ccResidual_ += asImp_().computeFluxForCellCenter(element, fvGeometry, elemVolVars, faceVars, scvf, elemFluxVarsCache[scvf]);
        }
    }

     /*!
     * \brief Evaluate the flux terms for face dofs
     */
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
                faceResiduals_[scvf.localFaceIdx()] += asImp_().computeFluxForFace(scvf, fvGeometry, elemVolVars, globalFaceVars);
        }
    }

     /*!
     * \brief Evaluate boundary conditions
     */
    void evalBoundary_(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const GlobalFaceVars& faceVars,
                       const ElementBoundaryTypes& bcTypes,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        DUNE_THROW(Dune::InvalidStateException,
           "The localResidual class does not provide "
           "a evalBoundary_() method.");
    }

     /*!
     * \brief Evaluate the volume term for a cell center dof for a stationary problem
     */
    template<class P = Problem>
    typename std::enable_if<Dumux::Capabilities::isStationary<P>::value, void>::type
    evalCCvolumeTerm_(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const SubControlVolume& scv,
                      const ElementVolumeVariables& prevElemVolVars,
                      const ElementVolumeVariables& curElemVolVars,
                      const GlobalFaceVars& prevFaceVars,
                      const GlobalFaceVars& curFaceVars,
                      const ElementBoundaryTypes &bcTypes)
    {
        const auto curExtrusionFactor = curElemVolVars[scv].extrusionFactor();

        // subtract the source term from the local rate
        CellCenterPrimaryVariables source = asImp_().computeSourceForCellCenter(element, fvGeometry, curElemVolVars, curFaceVars, scv);
        source *= scv.volume()*curExtrusionFactor;

        ccResidual_ -= source;
    }

     /*!
     * \brief Evaluate the volume term for a face dof for a stationary problem
     */
    template<class P = Problem>
    typename std::enable_if<Dumux::Capabilities::isStationary<P>::value, void>::type
    evalFaceVolumeTerm_(const Element &element,
                        const FVElementGeometry& fvGeometry,
                        const SubControlVolumeFace& scvf,
                        const ElementVolumeVariables& prevElemVolVars,
                        const ElementVolumeVariables& curElemVolVars,
                        const GlobalFaceVars& prevFaceVars,
                        const GlobalFaceVars& curFaceVars,
                        const ElementBoundaryTypes &bcTypes)
    {
        // the source term:
        auto faceSource = asImp_().computeSourceForFace(scvf, curElemVolVars, curFaceVars);
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        const auto curExtrusionFactor = curElemVolVars[scv].extrusionFactor();
        faceSource *= 0.5*scv.volume()*curExtrusionFactor;
        faceResiduals_[scvf.localFaceIdx()] -= faceSource;
    }

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
            evalCCvolumeTerm_(element, fvGeometry, scv, prevElemVolVars, curElemVolVars, prevFaceVars, curFaceVars, bcTypes);

            // now, treat the dofs on the facets:
            for(auto&& scvf : scvfs(fvGeometry))
            {
                evalFaceVolumeTerm_(element, fvGeometry, scvf, prevElemVolVars, curElemVolVars, prevFaceVars, curFaceVars, bcTypes);
            }
        }
    }

     /*!
     * \brief Evaluate the volume term for a cell center dof for a transient problem
     */
    template<class P = Problem>
    typename std::enable_if<!Dumux::Capabilities::isStationary<P>::value, void>::type
    evalCCvolumeTerm_(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const SubControlVolume& scv,
                      const ElementVolumeVariables& prevElemVolVars,
                      const ElementVolumeVariables& curElemVolVars,
                      const GlobalFaceVars& prevFaceVars,
                      const GlobalFaceVars& curFaceVars,
                      const ElementBoundaryTypes &bcTypes)
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

        ccStorageTerm_ = std::move(curCCStorage);
        ccStorageTerm_ -= std::move(prevCCStorage);
        ccStorageTerm_ *= scv.volume();
        ccStorageTerm_ /= problem().timeManager().timeStepSize();

        // add the storage term to the residual
        ccResidual_ += ccStorageTerm_;

        // subtract the source term from the local rate
        CellCenterPrimaryVariables source = asImp_().computeSourceForCellCenter(element, fvGeometry, curElemVolVars, curFaceVars, scv);
        source *= scv.volume()*curVolVars.extrusionFactor();

        ccResidual_ -= source;
    }

     /*!
     * \brief Evaluate the volume term for a face dof for a transient problem
     */
    template<class P = Problem>
    typename std::enable_if<!Dumux::Capabilities::isStationary<P>::value, void>::type
    evalFaceVolumeTerm_(const Element &element,
                        const FVElementGeometry& fvGeometry,
                        const SubControlVolumeFace& scvf,
                        const ElementVolumeVariables& prevElemVolVars,
                        const ElementVolumeVariables& curElemVolVars,
                        const GlobalFaceVars& prevFaceVars,
                        const GlobalFaceVars& curFaceVars,
                        const ElementBoundaryTypes &bcTypes)
    {
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& curVolVars = curElemVolVars[scv];
        const auto& prevVolVars = prevElemVolVars[scv];
        auto prevFaceStorage = asImp_().computeStorageForFace(scvf, prevVolVars, prevFaceVars);
        auto curFaceStorage = asImp_().computeStorageForFace(scvf, curVolVars, prevFaceVars);

        // the storage term
        faceStorageTerms_[scvf.localFaceIdx()] = std::move(curFaceStorage);
        faceStorageTerms_[scvf.localFaceIdx()] -= std::move(prevFaceStorage);
        faceStorageTerms_[scvf.localFaceIdx()] *= (scv.volume()/2.0);
        faceStorageTerms_[scvf.localFaceIdx()] /= problem().timeManager().timeStepSize();
        faceResiduals_[scvf.localFaceIdx()] += faceStorageTerms_[scvf.localFaceIdx()];

        // the source term:
        auto faceSource = asImp_().computeSourceForFace(scvf, curElemVolVars, curFaceVars);
        faceSource *= 0.5*scv.volume()*curVolVars.extrusionFactor();
        faceResiduals_[scvf.localFaceIdx()] -= faceSource;
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
        for (auto&& scv : scvs(fvGeometry))
        {
            // evaluate the volume terms (storage + source terms)
            evalCCvolumeTerm_(element, fvGeometry, scv, prevElemVolVars, curElemVolVars, prevFaceVars, curFaceVars, bcTypes);

            // now, treat the dofs on the facets:
            for(auto&& scvf : scvfs(fvGeometry))
            {
                evalFaceVolumeTerm_(element, fvGeometry, scvf, prevElemVolVars, curElemVolVars, prevFaceVars, curFaceVars, bcTypes);
            }
        }
    }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    CellCenterPrimaryVariables ccResidual_;
    CellCenterPrimaryVariables ccStorageTerm_;
    FaceSolutionVector faceResiduals_;
    FaceSolutionVector faceStorageTerms_;

private:
    Problem* problemPtr_;

};

}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH
