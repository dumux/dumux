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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_IMPLICIT_FLUXVARIABLES_HH
#define DUMUX_IMPLICIT_FLUXVARIABLES_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(NumComponents);
}

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the flux variables
 *        specializations are provided for combinations of physical processes
 */
template<class TypeTag, bool enableAdvection, bool enableMolecularDiffusion, bool enableEnergyBalance>
class FluxVariables {};


// specialization for pure advective flow (e.g. 1p/2p/3p immiscible darcy flow)
template<class TypeTag>
class FluxVariables<TypeTag, true, false, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Stencil = std::vector<IndexType>;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);

    enum
    {
        isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox),
        enableFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableFluxVariablesCache),
        constantBC = GET_PROP_VALUE(TypeTag, ConstantBoundaryConditions)
    };

public:
    FluxVariables() : problemPtr_(nullptr), scvFacePtr_(nullptr), boundaryVolVars_(nullptr)
    {}

    void update(const Problem& problem,
                const Element& element,
                const SubControlVolumeFace &scvFace)
    {
        problemPtr_ = &problem;
        scvFacePtr_ = &scvFace;

        // boundary vol vars only need to be handled at boundaries for cc methods
        if (scvFace.boundary() && !isBox)
            setBoundaryVolumeVariables_(problem, element, scvFace);

        // update the stencil if needed
        if (!enableFluxVarsCache && stencil_.empty())
            stencil_ = stencil(problem, scvFace);
    }



    template<typename FunctionType>
    Scalar advectiveFlux(const int phaseIdx, const FunctionType upwindFunction)
    {
        Scalar flux = AdvectionType::flux(problem(), scvFace(), phaseIdx, boundaryVolVars_);

        const auto* insideVolVars = &problem().model().curVolVars(scvFace().insideScvIdx());
        const VolumeVariables* outsideVolVars;

        if (scvFace().boundary())
        {
            if (boundaryVolVars_ == nullptr)
                DUNE_THROW(Dune::InvalidStateException, "Trying to access invalid boundary volume variables.");
            outsideVolVars = boundaryVolVars_;
        }
        else
            outsideVolVars = &problem().model().curVolVars(scvFace().outsideScvIdx());

        if (std::signbit(flux))
            return flux*upwindFunction(*outsideVolVars, *insideVolVars);
        else
            return flux*upwindFunction(*insideVolVars, *outsideVolVars);

        return flux;
    }

    // interface allowing for stencil information without having to update the flux vars.
    // this becomes useful when the element is not at hand for a call to update(...), e.g. during the assembly - localjacobian.hh
    template <typename T = TypeTag>
    const typename std::enable_if<GET_PROP_VALUE(T, EnableFluxVariablesCache), Stencil>::type& stencil(const Problem& problem, const SubControlVolumeFace& scvFace) const
    { return problem.model().fluxVarsCache(scvFace).stencil(); }

    // provide interface in case caching is disabled
    template <typename T = TypeTag>
    const typename std::enable_if<!GET_PROP_VALUE(T, EnableFluxVariablesCache), Stencil>::type& stencil(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        if (stencil_.empty())
            stencil_ = computeFluxStencil(problem, scvFace);
        return stencil_;
    }

    // returns the boundary vol vars belonging to this flux variables object
    VolumeVariables getBoundaryVolumeVariables(const Problem& problem, const Element& element, const SubControlVolumeFace& scvFace)
    {
        VolumeVariables boundaryVolVars;

        const auto insideScvIdx = scvFace.insideScvIdx();
        const auto& insideScv = problem.model().fvGeometries().subControlVolume(insideScvIdx);
        const auto dirichletPriVars = problem.dirichlet(element, scvFace);
        boundaryVolVars.update(dirichletPriVars, problem, element, insideScv);

        return boundaryVolVars;
    }

    const Problem& problem() const
    {
        if (problemPtr_ == nullptr)
            DUNE_THROW(Dune::InvalidStateException, "Problem pointer not valid. Call update before using the flux variables class.");
        return *problemPtr_;
    }

    const SubControlVolumeFace& scvFace() const
    {
        if (scvFacePtr_ == nullptr)
            DUNE_THROW(Dune::InvalidStateException, "Scv face pointer not valid. Call update before using the flux variables class.");
        return *scvFacePtr_;
    }

    // when caching is enabled, get the stencil from the cache class
    template <typename T = TypeTag>
    const typename std::enable_if<GET_PROP_VALUE(T, EnableFluxVariablesCache), Stencil>::type& stencil() const
    {
        if (problemPtr_ == nullptr || scvFacePtr_ == nullptr)
            DUNE_THROW(Dune::InvalidStateException, "Calling the stencil() method before update of the FluxVariables object.");
        return problem().model().fluxVarsCache(scvFace()).stencil();
    }

    // when caching is disabled, return the private stencil variable. The update(...) routine has to be called beforehand.
    template <typename T = TypeTag>
    const typename std::enable_if<!GET_PROP_VALUE(T, EnableFluxVariablesCache), Stencil>::type& stencil()
    {
        if (stencil_.empty())
            DUNE_THROW(Dune::InvalidStateException, "Calling the stencil() method before update of the FluxVariables object.");
        return stencil_;
    }

    Stencil computeFluxStencil(const Problem& problem, const SubControlVolumeFace& scvFace)
    { return AdvectionType::stencil(problem, scvFace); }

private:

    // if flux variables caching is enabled, we use the boundary volume variables stored in the cache class
    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableFluxVariablesCache) && GET_PROP_VALUE(T, ConstantBoundaryConditions)>::type
    setBoundaryVolumeVariables_(const Problem& problem, const Element& element, const SubControlVolumeFace& scvFace)
    {
        boundaryVolVars_ = &problem.model().fluxVarsCache(scvFace).boundaryVolumeVariables();
    }

    // if flux variables caching is disabled, we update the boundary volume variables
    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, ConstantBoundaryConditions) || !GET_PROP_VALUE(T, EnableFluxVariablesCache)>::type
    setBoundaryVolumeVariables_(const Problem& problem, const Element& element, const SubControlVolumeFace& scvFace)
    {
        VolumeVariables tmp;

        const auto dirichletPriVars = problem.dirichlet(element, scvFace);
        const auto insideScvIdx = scvFace.insideScvIdx();
        const auto& insideScv = problem.model().fvGeometries().subControlVolume(insideScvIdx);
        tmp.update(dirichletPriVars, problem, element, insideScv);

        boundaryVolVars_ = new VolumeVariables(std::move(tmp));
    }

    const Problem *problemPtr_;              //! Pointer to the problem
    const SubControlVolumeFace *scvFacePtr_; //! Pointer to the sub control volume face for which the flux variables are created
    const VolumeVariables* boundaryVolVars_; //! boundary volume variables in case of Dirichlet boundaries
    // the flux stencil
    Stencil stencil_;
};


// specialization for isothermal advection molecularDiffusion equations
template<class TypeTag>
class FluxVariables<TypeTag, true, true, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Stencil = std::vector<IndexType>;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:

    FluxVariables()
    {}

    void update(const Problem& problem,
                const Element& element,
                const SubControlVolumeFace &scvFace)
    {
        problemPtr_ = &problem;
        scvFacePtr_ = &scvFace;

        if (scvFace.boundary())
        {
            if(!boundaryVolVars_)
                boundaryVolVars_ = Dune::Std::make_unique<VolumeVariables>();
        }

        boundaryVolVarsUpdated_ = false;

        // TODO compute stencil as unity of advective and diffusive flux stencils
        stencil_ = AdvectionType::stencil(scvFace);
    }

    void beginFluxComputation()
    {
        // TODO if constant BC, do not set to false
        boundaryVolVarsUpdated_ = false;

        // calculate advective fluxes and store them
        advectiveVolFluxes_ = new std::array<Scalar, numPhases>();
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            (*advectiveVolFluxes_)[phaseIdx] = AdvectionType::flux(problem(), scvFace(), phaseIdx, boundaryVolVars_.get(), boundaryVolVarsUpdated_);
            boundaryVolVarsUpdated_ = true;
        }
    }

    void endFluxComputation()
    {
        delete advectiveVolFluxes_;
    }

    template<typename FunctionType>
    Scalar advectiveFlux(const int phaseIdx, const FunctionType upwindFunction)
    {
        const auto* insideVolVars = &problem().model().curVolVars(scvFace().insideScvIdx());
        const VolumeVariables* outsideVolVars;

        if (scvFace().boundary())
            outsideVolVars = boundaryVolVars_.get();
        else
            outsideVolVars = &problem().model().curVolVars(scvFace().outsideScvIdx());

        if (std::signbit((*advectiveVolFluxes_)[phaseIdx]))
            return (*advectiveVolFluxes_)[phaseIdx]*upwindFunction(*outsideVolVars, *insideVolVars);
        else
            return (*advectiveVolFluxes_)[phaseIdx]*upwindFunction(*insideVolVars, *outsideVolVars);
    }

    Scalar molecularDiffusionFlux(const int phaseIdx, const int compIdx)
    {
        Scalar flux = MolecularDiffusionType::flux(problem(), scvFace(), phaseIdx, compIdx, boundaryVolVars_.get(), boundaryVolVarsUpdated_);
        boundaryVolVarsUpdated_ = true;
        return flux;
    }

    const Stencil& stencil() const
    {
        return stencil_;
    }

    const Problem& problem() const
    {
        return *problemPtr_;
    }

    const SubControlVolumeFace& scvFace() const
    {
        return *scvFacePtr_;
    }

private:
    const Problem *problemPtr_;              //! Pointer to the problem
    const SubControlVolumeFace *scvFacePtr_; //! Pointer to the sub control volume face for which the flux variables are created

    // storage for calculated advective fluxes to not having to calculate them again
    std::array<Scalar, numPhases>* advectiveVolFluxes_;
    // boundary volume variables in case of Dirichlet boundaries
    std::unique_ptr<VolumeVariables> boundaryVolVars_;
    bool boundaryVolVarsUpdated_;
    // the flux stencil
    Stencil stencil_;
};


// specialization for pure molecularDiffusion_
template<class TypeTag>
class FluxVariables<TypeTag, false, true, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::set<IndexType>;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scvf)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                molecularDiffusion_.update(problem, scvf, phaseIdx, compIdx);
    }

    const MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx) const
    {
        return molecularDiffusion_[phaseIdx][compIdx];
    }

    MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx)
    {
        return molecularDiffusion_[phaseIdx][compIdx];
    }

private:
    std::array< std::array<MolecularDiffusionType, numComponents>, numPhases> molecularDiffusion_;
};


// specialization for non-isothermal advective flow (e.g. non-isothermal one-phase darcy equation)
template<class TypeTag>
class FluxVariables<TypeTag, true, false, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::set<IndexType>;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

public:
    void update(const Problem& problem, const SubControlVolumeFace &scvf)
    {
        advection_.update(problem, scvf);
        heatConduction_.update(problem, scvf);
    }

    const AdvectionType& advection() const
    {
        return advection_;
    }

    AdvectionType& advection()
    {
        return advection_;
    }

    const HeatConductionType& heatConduction() const
    {
        return heatConduction_;
    }

    HeatConductionType& heatConduction()
    {
        return heatConduction_;
    }

private:
    AdvectionType advection_;
    HeatConductionType heatConduction_;
};


// specialization for non-isothermal advection molecularDiffusion_ equations
template<class TypeTag>
class FluxVariables<TypeTag, true, true, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::set<IndexType>;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scvf)
    {
        advection_.update(problem, scvf);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                molecularDiffusion_.update(problem, scvf, phaseIdx, compIdx);
        heatConduction_.update(problem, scvf);
    }

    const AdvectionType& advection() const
    {
        return advection_;
    }

    AdvectionType& advection()
    {
        return advection_;
    }

    const HeatConductionType& heatConduction() const
    {
        return heatConduction_;
    }

    HeatConductionType& heatConduction()
    {
        return heatConduction_;
    }

    const MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx) const
    {
        return molecularDiffusion_[phaseIdx][compIdx];
    }

    MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx)
    {
        return molecularDiffusion_[phaseIdx][compIdx];
    }

private:
    AdvectionType advection_;
    std::array< std::array<MolecularDiffusionType, numComponents>, numPhases> molecularDiffusion_;
    HeatConductionType heatConduction_;
};


// specialization for non-isothermal molecularDiffusion_
template<class TypeTag>
class FluxVariables<TypeTag, false, true, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::set<IndexType>;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    void update(const Problem& problem, const SubControlVolumeFace &scvf)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                molecularDiffusion_.update(problem, scvf, phaseIdx, compIdx);
        heatConduction_.update(problem, scvf);
    }

    const HeatConductionType& heatConduction() const
    {
        return heatConduction_;
    }

    HeatConductionType& heatConduction()
    {
        return heatConduction_;
    }

    const MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx) const
    {
        return molecularDiffusion_[phaseIdx][compIdx];
    }

    MolecularDiffusionType& molecularDiffusion(const int phaseIdx, const int compIdx)
    {
        return molecularDiffusion_[phaseIdx][compIdx];
    }

private:
    std::array< std::array<MolecularDiffusionType, numComponents>, numPhases> molecularDiffusion_;
    HeatConductionType heatConduction_;
};


// specialization for pure heat conduction (e.g. the heat equation)
template<class TypeTag>
class FluxVariables<TypeTag, false, false, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::set<IndexType>;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

public:
    void update(const Problem& problem, const SubControlVolumeFace &scvf)
    {
        heatConduction_.update(problem, scvf);
    }

    const HeatConductionType& heatConduction() const
    {
        return heatConduction_;
    }

    HeatConductionType& heatConduction()
    {
        return heatConduction_;
    }

private:
    HeatConductionType heatConduction_;
};

} // end namespace

#endif
