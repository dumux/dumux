// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredElementVolumeVariables
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_ELEMENT_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_ELEMENT_VOLUMEVARIABLES_HH

#include <algorithm>
#include <cassert>
#include <vector>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/common/typetraits/vector.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the element volume variables vector for the staggered model
 */
template<class GVV, bool cachingEnabled>
class StaggeredElementVolumeVariables
{};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for the element volume variables vector for the staggered model.
          Specialization in case the volume variables are stored globally.
 */
template<class GVV>
class StaggeredElementVolumeVariables<GVV, /*cachingEnabled*/true>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    StaggeredElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars)
    , numScv_(gridVolVars.problem().gridGeometry().numScv())
    {}

    //! operator for the access with an scv
    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    {
        if (scv.dofIndex() < numScv_)
            return gridVolVars().volVars(scv.dofIndex());
        else
            return boundaryVolumeVariables_[getLocalIdx_(scv.dofIndex())];
    }

    //! operator for the access with an index
    //! needed for Staggered methods for the access to the boundary volume variables
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
    template<class FVElementGeometry, class SolutionVector>
    StaggeredElementVolumeVariables bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                         const FVElementGeometry& fvGeometry,
                                         const SolutionVector& sol) &&
    {
        this->bind_(element, fvGeometry, sol);
        return std::move(*this);
    }

    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol) &
    { this->bind_(element, fvGeometry, sol); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    StaggeredElementVolumeVariables bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                                const FVElementGeometry& fvGeometry,
                                                const SolutionVector& sol) &&
    {
        this->bindElement_(element, fvGeometry, sol);
        return std::move(*this);
    }

    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol) &
    { this->bindElement_(element, fvGeometry, sol); }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:

    //! Clear all local storage
    void clear_()
    {
        boundaryVolVarIndices_.clear();
        boundaryVolumeVariables_.clear();
    }

    //! Binding of an element, prepares the volume variables within the element stencil
    //! called by the local jacobian to prepare element assembly. Specialization callable with MultiTypeBlockVector.
    template<class FVElementGeometry, class SolutionVector, typename std::enable_if_t<isMultiTypeBlockVector<SolutionVector>::value, int> = 0>
    void bind_(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
               const FVElementGeometry& fvGeometry,
               const SolutionVector& sol)
    {
        // forward to the actual method
        bind_(element, fvGeometry, sol[FVElementGeometry::GridGeometry::cellCenterIdx()]);
    }

    //! Binding of an element, prepares the volume variables within the element stencil
    //! called by the local jacobian to prepare element assembly
    template<class FVElementGeometry, class SolutionVector, typename std::enable_if_t<!isMultiTypeBlockVector<SolutionVector>::value, int> = 0>
    void bind_(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
               const FVElementGeometry& fvGeometry,
               const SolutionVector& sol)
    {
        if (!fvGeometry.hasBoundaryScvf())
            return;

        clear_();
        boundaryVolVarIndices_.reserve(fvGeometry.numScvf());
        boundaryVolumeVariables_.reserve(fvGeometry.numScvf());

        // handle the boundary volume variables
        for (auto&& scvf : scvfs(fvGeometry))
        {
            // if we are not on a boundary, skip the rest
            if (!scvf.boundary())
                continue;

            const auto& problem = gridVolVars().problem();
            auto boundaryPriVars = gridVolVars().getBoundaryPriVars(problem, sol, element, scvf);
            const auto elemSol = elementSolution<FVElementGeometry>(std::move(boundaryPriVars));
            auto&& scvI = fvGeometry.scv(scvf.insideScvIdx());

            VolumeVariables volVars;
            volVars.update(elemSol,
                           problem,
                           element,
                           scvI);

           boundaryVolumeVariables_.emplace_back(std::move(volVars));
           boundaryVolVarIndices_.push_back(scvf.outsideScvIdx());
        }
    }

    //! function to prepare the vol vars within the element
    template<class FVElementGeometry, class SolutionVector>
    void bindElement_(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                      const FVElementGeometry& fvGeometry,
                      const SolutionVector& sol)
    {}

    const GridVolumeVariables* gridVolVarsPtr_;

    //! map a global scv index to the local storage index
    int getLocalIdx_(const int volVarIdx) const
    {
        auto it = std::find(boundaryVolVarIndices_.begin(), boundaryVolVarIndices_.end(), volVarIdx);
        assert(it != boundaryVolVarIndices_.end() && "Could not find the current volume variables for volVarIdx!");
        return std::distance(boundaryVolVarIndices_.begin(), it);
    }

    std::vector<std::size_t> boundaryVolVarIndices_;
    std::vector<VolumeVariables> boundaryVolumeVariables_;
    const std::size_t numScv_;
};


/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for the element volume variables vector for the staggered model.
          Specialization in case the volume variables are not stored globally.
 */
template<class GVV>
class StaggeredElementVolumeVariables<GVV, /*cachingEnabled*/false>
{
    using PrimaryVariables = typename GVV::VolumeVariables::PrimaryVariables;

public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    StaggeredElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    StaggeredElementVolumeVariables bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                         const FVElementGeometry& fvGeometry,
                                         const SolutionVector& sol) &&
    {
        this->bind_(element, fvGeometry, sol);
        return std::move(*this);
    }

    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol) &
    { this->bind_(element, fvGeometry, sol); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    StaggeredElementVolumeVariables bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                                const FVElementGeometry& fvGeometry,
                                                const SolutionVector& sol) &&
    {
        this->bindElement_(element, fvGeometry, sol);
        return std::move(*this);
    }

    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol) &
    { this->bindElement_(element, fvGeometry, sol); }

    //! const operator for the access with an scv
    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return volumeVariables_[getLocalIdx_(scv.dofIndex())]; }

    //! operator for the access with an scv
    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    VolumeVariables& operator [](const SubControlVolume& scv)
    { return volumeVariables_[getLocalIdx_(scv.dofIndex())]; }

    //! const operator for the access with an index
    const VolumeVariables& operator [](std::size_t scvIdx) const
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    //! operator for the access with an index
    VolumeVariables& operator [](std::size_t scvIdx)
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
        //! Binding of an element, prepares the volume variables within the element stencil
    //! called by the local jacobian to prepare element assembly. Specialization callable with MultiTypeBlockVector.
    template<class FVElementGeometry, class SolutionVector, typename std::enable_if_t<isMultiTypeBlockVector<SolutionVector>::value, int> = 0>
    void bind_(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
               const FVElementGeometry& fvGeometry,
               const SolutionVector& sol)
    {
        // forward to the actual method
        bind_(element, fvGeometry, sol[FVElementGeometry::GridGeometry::cellCenterIdx()]);
    }

    //! Binding of an element, prepares the volume variables within the element stencil
    //! called by the local jacobian to prepare element assembly
    template<class FVElementGeometry, class SolutionVector, typename std::enable_if_t<!isMultiTypeBlockVector<SolutionVector>::value, int> = 0>
    void bind_(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
               const FVElementGeometry& fvGeometry,
               const SolutionVector& sol)
    {
        clear_();

        const auto& problem = gridVolVars().problem();
        const auto& gridGeometry = fvGeometry.gridGeometry();
        const auto globalI = gridGeometry.elementMapper().index(element);
        const auto& map = gridGeometry.connectivityMap();
        constexpr auto cellCenterIdx = FVElementGeometry::GridGeometry::cellCenterIdx();
        const auto& connectivityMapI = map(cellCenterIdx, cellCenterIdx, globalI);
        const auto numDofs = connectivityMapI.size();

        auto&& scvI = fvGeometry.scv(globalI);

        // resize local containers to the required size (for internal elements)
        volumeVariables_.resize(numDofs+1);
        volVarIndices_.resize(numDofs+1);
        int localIdx = 0;

        // Lambda to update the volume variables of the given index
        auto doVolVarUpdate = [&](int globalJ)
        {
            const auto& elementJ = gridGeometry.element(globalJ);
            auto&& scvJ = fvGeometry.scv(globalJ);
            const auto elemSol = makeElementSolutionFromCellCenterPrivars<PrimaryVariables>(sol[globalJ]);
            volumeVariables_[localIdx].update(elemSol,
                                              problem,
                                              elementJ,
                                              scvJ);
            volVarIndices_[localIdx] = scvJ.dofIndex();
            ++localIdx;
        };

        // Update the volume variables of the element at hand
        doVolVarUpdate(globalI);

        // Update the volume variables of the neighboring elements
        for (const auto& globalJ : connectivityMapI)
            doVolVarUpdate(globalJ);

        if (fvGeometry.hasBoundaryScvf())
        {
            // Update boundary volume variables
            for (auto&& scvf : scvfs(fvGeometry))
            {
                // if we are not on a boundary, skip to the next scvf
                if (!scvf.boundary())
                    continue;

                volumeVariables_.resize(localIdx+1);
                volVarIndices_.resize(localIdx+1);

                auto boundaryPriVars = gridVolVars().getBoundaryPriVars(problem, sol, element, scvf);
                auto elemSol = elementSolution<FVElementGeometry>(std::move(boundaryPriVars));
                volumeVariables_[localIdx].update(elemSol,
                                                  problem,
                                                  element,
                                                  scvI);
                volVarIndices_[localIdx] = scvf.outsideScvIdx();
                ++localIdx;
            }
        }
    }

    //! Binding of an element, prepares the volume variables within the element stencil
    //! called by the local jacobian to prepare element assembly. Specialization callable with MultiTypeBlockVector.
    template<class FVElementGeometry, class SolutionVector, typename std::enable_if_t<isMultiTypeBlockVector<SolutionVector>::value, int> = 0>
    void bindElement_(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                      const FVElementGeometry& fvGeometry,
                      const SolutionVector& sol)
    {
        // forward to the actual method
        bindElement_(element, fvGeometry, sol[FVElementGeometry::GridGeometry::cellCenterIdx()]);
    }

    //! Binding of an element, prepares only the volume variables of the element.
    //! Specialization for Staggered models
    template<class FVElementGeometry, class SolutionVector, typename std::enable_if_t<!isMultiTypeBlockVector<SolutionVector>::value, int> = 0>
    void bindElement_(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                      const FVElementGeometry& fvGeometry,
                      const SolutionVector& sol)
    {
        clear_();

        const auto globalI = fvGeometry.gridGeometry().elementMapper().index(element);
        volumeVariables_.resize(1);
        volVarIndices_.resize(1);

        // update the volume variables of the element
        auto&& scv = fvGeometry.scv(globalI);

        const auto elemSol = makeElementSolutionFromCellCenterPrivars<PrimaryVariables>(sol[globalI]);
        volumeVariables_[0].update(elemSol,
                                   gridVolVars().problem(),
                                   element,
                                   scv);
        volVarIndices_[0] = scv.dofIndex();
    }

    //! Clear all local storage
    void clear_()
    {
        volVarIndices_.clear();
        volumeVariables_.clear();
    }

    const GridVolumeVariables* gridVolVarsPtr_;

    int getLocalIdx_(const int volVarIdx) const
    {
        auto it = std::find(volVarIndices_.begin(), volVarIndices_.end(), volVarIdx);
        assert(it != volVarIndices_.end() && "Could not find the current volume variables for volVarIdx!");
        return std::distance(volVarIndices_.begin(), it);
    }

    std::vector<std::size_t> volVarIndices_;
    std::vector<VolumeVariables> volumeVariables_;
};

} // end namespace Dumux

#endif
