// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredElementVolumeVariables
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_ELEMENTVOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_ELEMENTVOLUMEVARIABLES_HH

#include <algorithm>
#include <cassert>
#include <vector>
#include <utility>

#include <dumux/discretization/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Base class for the face variables vector
 */
template<class GridVolumeVariables, bool cachingEnabled>
class FaceCenteredStaggeredElementVolumeVariables;

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Class for the face variables vector. Specialization for the case of storing the face variables globally.
 */
template<class GVV>
class FaceCenteredStaggeredElementVolumeVariables<GVV, /*cachingEnabled*/true>
{
    using ThisType = FaceCenteredStaggeredElementVolumeVariables<GVV, /*cachingEnabled*/true>;
    using GridGeometry = std::decay_t<decltype(std::declval<GVV>().problem().gridGeometry())>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    FaceCenteredStaggeredElementVolumeVariables(const GridVolumeVariables& gridVolumeVariables)
    : gridVolumeVariablesPtr_(&gridVolumeVariables)
    , numScv_(gridVolumeVariables.problem().gridGeometry().numScv())
    {}

    //! operator for the access with an scvf
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    {
        if (scv.index() < numScv_)
            return gridVolVars().volVars(scv.index());
        else
            return boundaryVolumeVariables_[getLocalIdx_(scv.index())];
    }

    //! operator for the access with an index
    //! needed for cc methods for the access to the boundary volume variables
    const VolumeVariables& operator [](const std::size_t scvIdx) const
    {
        if (scvIdx < numScv_)
            return gridVolVars().volVars(scvIdx);
        else
            return boundaryVolumeVariables_[getLocalIdx_(scvIdx)];
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class SolutionVector>
    ThisType bind(const typename FVElementGeometry::Element& element,
                  const FVElementGeometry& fvGeometry,
                  const SolutionVector& sol) &&
    {
        this->bind(element, fvGeometry, sol);
        return std::move(*this);
    }

    //! For compatibility reasons with the case of not storing the face vars.
    //! function to be called before assembling an element, preparing the vol vars within the stencil
    template<class SolutionVector>
    void bind(const typename FVElementGeometry::Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol) &
    {
        if (!fvGeometry.hasBoundaryScvf())
            return;

        clear_();
        // upper bound of size is the number of all scvfs minus frontal scvfs
        boundaryVolVarIndices_.reserve(fvGeometry.numScvf()-element.subEntities(1));
        boundaryVolumeVariables_.reserve(fvGeometry.numScvf()-element.subEntities(1));

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary() || scvf.isFrontal() || scvf.processorBoundary())
                continue;

            // check if boundary is a pure dirichlet boundary
            const auto& problem = gridVolVars().problem();
            const auto bcTypes = problem.boundaryTypes(element, scvf);

            auto addBoundaryVolVars = [&](const auto& scvFace)
            {
                const auto& scvI = fvGeometry.scv(scvFace.insideScvIdx());
                typename VolumeVariables::PrimaryVariables pv(
                    problem.dirichlet(element, scvFace)[scvI.dofAxis()]
                );
                const auto dirichletPriVars = elementSolution<FVElementGeometry>(pv);

                VolumeVariables volVars;
                volVars.update(dirichletPriVars, problem, element, scvI);

                boundaryVolumeVariables_.emplace_back(std::move(volVars));
                boundaryVolVarIndices_.push_back(scvFace.outsideScvIdx());
            };

            if (bcTypes.hasDirichlet())
            {
                addBoundaryVolVars(scvf);
                continue;
            }

            // treat domain corners
            if (const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf); orthogonalScvf.boundary())
                if (problem.boundaryTypes(element, orthogonalScvf).hasDirichlet())
                    addBoundaryVolVars(scvf);

        }

        assert(boundaryVolumeVariables_.size() == boundaryVolVarIndices_.size());
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class SolutionVector>
    ThisType bindElement(const typename FVElementGeometry::Element& element,
                         const FVElementGeometry& fvGeometry,
                         const SolutionVector& sol) &&
    {
        this->bindElement(element, fvGeometry, sol);
        return std::move(*this);
    }
    //! Binding of an element, prepares only the face variables of the element
    //! specialization for Staggered models
    template<class SolutionVector>
    void bindElement(const typename FVElementGeometry::Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol) &
    {}


    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolumeVariablesPtr_; }

    //! Returns true if volVars exist for the given scv index
    bool hasVolVars(const std::size_t scvIdx) const
    {
        if (scvIdx < numScv_)
            return true;
        else
        {
            const auto it = std::find(boundaryVolVarIndices_.begin(), boundaryVolVarIndices_.end(), scvIdx);
            return it != boundaryVolVarIndices_.end();
        }
    }

private:
    //! Clear all local storage
    void clear_()
    {
        boundaryVolVarIndices_.clear();
        boundaryVolumeVariables_.clear();
    }

    //! map a global scv index to the local storage index
    int getLocalIdx_(const std::size_t volVarIdx) const
    {
        const auto it = std::find(boundaryVolVarIndices_.begin(), boundaryVolVarIndices_.end(), volVarIdx);
        assert(it != boundaryVolVarIndices_.end() && "Could not find the current volume variables for volVarIdx!");
        return std::distance(boundaryVolVarIndices_.begin(), it);
    }

    const GridVolumeVariables* gridVolumeVariablesPtr_;
    const std::size_t numScv_;
    std::vector<std::size_t> boundaryVolVarIndices_;
    std::vector<VolumeVariables> boundaryVolumeVariables_;
};

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Class for the face variables vector. Specialization for the case of not storing the face variables globally.
 */
template<class GVV>
class FaceCenteredStaggeredElementVolumeVariables<GVV, /*cachingEnabled*/false>
{
    using ThisType = FaceCenteredStaggeredElementVolumeVariables<GVV, /*cachingEnabled*/false>;
    using GridGeometry = std::decay_t<decltype(std::declval<GVV>().problem().gridGeometry())>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

    static constexpr auto dim = GridGeometry::GridView::dimension;
    static constexpr auto numInsideVolVars = dim * 2;
    static constexpr auto numOutsideVolVars = numInsideVolVars * 2 * (dim - 1);

public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    FaceCenteredStaggeredElementVolumeVariables(const GridVolumeVariables& globalFacesVars)
    : gridVolumeVariablesPtr_(&globalFacesVars) {}

    //! const operator for the access with an scvf
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return volumeVariables_[getLocalIdx_(scv.index())]; }

    //! const operator for the access with an index
    const VolumeVariables& operator [](const std::size_t scvIdx) const
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    //! operator for the access with an scvf
    VolumeVariables& operator [](const SubControlVolume& scv)
    { return volumeVariables_[getLocalIdx_(scv.index())]; }

    // operator for the access with an index
    VolumeVariables& operator [](const std::size_t scvIdx)
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class SolutionVector>
    ThisType bind(const typename FVElementGeometry::Element& element,
                  const FVElementGeometry& fvGeometry,
                  const SolutionVector& sol) &&
    {
        this->bind_(element, fvGeometry, sol);
        return std::move(*this);
    }

    template<class SolutionVector>
    void bind(const typename FVElementGeometry::Element& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol) &
    {
        this->bind_(element, fvGeometry, sol);
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class SolutionVector>
    ThisType bindElement(const typename FVElementGeometry::Element& element,
                         const FVElementGeometry& fvGeometry,
                         const SolutionVector& sol) &&
    {
        this->bindElement_(element, fvGeometry, sol);
        return std::move(*this);
    }

    template<class SolutionVector>
    void bindElement(const typename FVElementGeometry::Element& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol) &
    { this->bindElement_(element, fvGeometry, sol); }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolumeVariablesPtr_; }

    //! Returns true if volVars exist for the given scv index
    bool hasVolVars(const std::size_t scvIdx) const
    { return volVarsInserted_(scvIdx); }

private:
    //! For compatibility reasons with the case of not storing the vol vars.
    //! function to be called before assembling an element, preparing the vol vars within the stencil
    template<class SolutionVector>
    void bind_(const typename FVElementGeometry::Element& element,
               const FVElementGeometry& fvGeometry,
               const SolutionVector& sol)
    {
        clear_();

        const auto& problem = gridVolVars().problem();
        const auto& gridGeometry = fvGeometry.gridGeometry();

        volVarIndices_.reserve(numInsideVolVars + numInsideVolVars);
        volumeVariables_.reserve(numInsideVolVars + numInsideVolVars);

        for (const auto& scv : scvs(fvGeometry))
        {
            for (const auto otherScvIdx : gridGeometry.connectivityMap()[scv.index()])
            {
                if (!volVarsInserted_(otherScvIdx))
                {
                    const auto& otherScv = fvGeometry.scv(otherScvIdx);
                    volVarIndices_.push_back(otherScvIdx);
                    volumeVariables_.emplace_back();
                    const auto& otherElement = gridGeometry.element(otherScv.elementIndex());
                    volumeVariables_.back().update(
                        elementSolution(otherElement, sol, gridGeometry),
                        problem, otherElement, otherScv
                    );
                }
            }
        }

        if (fvGeometry.hasBoundaryScvf())
        {
            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (!scvf.boundary() || scvf.isFrontal())
                    continue;

                // check if boundary is a pure dirichlet boundary
                const auto& problem = gridVolVars().problem();
                const auto bcTypes = problem.boundaryTypes(element, scvf);

                auto addBoundaryVolVars = [&](const auto& scvFace)
                {
                    const auto& scvI = fvGeometry.scv(scvFace.insideScvIdx());
                    typename VolumeVariables::PrimaryVariables pv(
                        problem.dirichlet(element, scvFace)[scvI.dofAxis()]
                    );
                    const auto dirichletPriVars = elementSolution<FVElementGeometry>(pv);

                    VolumeVariables volVars;
                    volVars.update(dirichletPriVars,
                                   problem,
                                   element,
                                   scvI);

                    volumeVariables_.emplace_back(std::move(volVars));
                    volVarIndices_.push_back(scvFace.outsideScvIdx());
                };

                if (bcTypes.hasDirichlet())
                {
                    addBoundaryVolVars(scvf);
                    continue;
                }

                // treat domain corners
                if (const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf); orthogonalScvf.boundary())
                    if (problem.boundaryTypes(element, orthogonalScvf).hasDirichlet())
                        addBoundaryVolVars(scvf);

            }
        }
    }

    //! Binding of an element, prepares only the face variables of the element
    //! specialization for Staggered models
    template<class SolutionVector>
    void bindElement_(const typename FVElementGeometry::Element& element,
                      const FVElementGeometry& fvGeometry,
                      const SolutionVector& sol)
    {
        clear_();
        const auto& problem = gridVolVars().problem();
        const auto& gridGeometry = fvGeometry.gridGeometry();
        volVarIndices_.reserve(numInsideVolVars);

        for (const auto& scv : scvs(fvGeometry))
        {
            volVarIndices_.push_back(scv.index());
            volumeVariables_.emplace_back();
            volumeVariables_.back().update(
                elementSolution(element, sol, gridGeometry),
                problem, element, scv
            );
        }
    }

    //! Clear all local storage
    void clear_()
    {
        volVarIndices_.clear();
        volumeVariables_.clear();
    }

    bool volVarsInserted_(const std::size_t scvIdx) const
    {
        return std::find(volVarIndices_.begin(), volVarIndices_.end(), scvIdx) != volVarIndices_.end();
    }

    int getLocalIdx_(const int scvfIdx) const
    {
        const auto it = std::find(volVarIndices_.begin(), volVarIndices_.end(), scvfIdx);
        assert(it != volVarIndices_.end() && "Could not find the current face variables for scvfIdx!");
        return std::distance(volVarIndices_.begin(), it);
    }

    const GridVolumeVariables* gridVolumeVariablesPtr_;
    std::vector<std::size_t> volVarIndices_;
    std::vector<VolumeVariables> volumeVariables_;
};

} // end namespace Dumux

#endif
