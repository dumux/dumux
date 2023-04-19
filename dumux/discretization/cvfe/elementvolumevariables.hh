// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief The local volume variables class
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_ELEMENT_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_CVFE_ELEMENT_VOLUMEVARIABLES_HH

#include <type_traits>
#include <utility>
#include <vector>

#include <dumux/discretization/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup CVFEDiscretization
 * \brief The local (stencil) volume variables class for control-volume finite element
 * \note The class is specialized for versions with and without caching
 * \tparam GVV the grid volume variables type
 * \tparam cachingEnabled if the cache is enabled
 */
template<class GVV, bool cachingEnabled>
class CVFEElementVolumeVariables;

/*!
 * \ingroup CVFEDiscretization
 * \brief The local (stencil) volume variables class for control-volume finite element with caching
 * \note the volume variables are stored for the whole grid view in the corresponding GridVolumeVariables class
 */
template<class GVV>
class CVFEElementVolumeVariables<GVV, /*cachingEnabled*/true>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    CVFEElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    const VolumeVariables& operator [](std::size_t scvIdx) const
    { return gridVolVars().volVars(eIdx_, scvIdx); }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return gridVolVars().volVars(eIdx_, scv.indexInElement()); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    CVFEElementVolumeVariables bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                   const FVElementGeometry& fvGeometry,
                                   const SolutionVector& sol)  &&
    {
        this->bindElement(element, fvGeometry, sol);
        return std::move(*this);
    }

    // For compatibility reasons with the case of not storing the vol vars.
    // function to be called before assembling an element, preparing the vol vars within the stencil
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol) &
    {
        bindElement(element, fvGeometry, sol);
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    CVFEElementVolumeVariables bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                          const FVElementGeometry& fvGeometry,
                                          const SolutionVector& sol)  &&
    {
        this->bindElement(element, fvGeometry, sol);
        return std::move(*this);
    }

    // function to prepare the vol vars within the element
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol) &
    {
        eIdx_ = fvGeometry.gridGeometry().elementMapper().index(element);
    }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    const GridVolumeVariables* gridVolVarsPtr_;
    std::size_t eIdx_;
};


/*!
 * \ingroup CVFEDiscretization
 * \brief The local (stencil) volume variables class for control-volume finite element without caching
 */
template<class GVV>
class CVFEElementVolumeVariables<GVV, /*cachingEnabled*/false>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! Constructor
    CVFEElementVolumeVariables(const GridVolumeVariables& gridVolVars)
    : gridVolVarsPtr_(&gridVolVars) {}

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    CVFEElementVolumeVariables bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                   const FVElementGeometry& fvGeometry,
                                   const SolutionVector& sol)  &&
    {
        this->bindElement(element, fvGeometry, sol);
        return std::move(*this);
    }

    // specialization for control-volume finite element, simply forwards to the bindElement method
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol) &
    {
        bindElement(element, fvGeometry, sol);
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    CVFEElementVolumeVariables bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                                          const FVElementGeometry& fvGeometry,
                                          const SolutionVector& sol)  &&
    {
        this->bindElement(element, fvGeometry, sol);
        return std::move(*this);
    }

    // specialization for control-volume finite element
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol) &
    {
        // get the solution at the dofs of the element
        auto elemSol = elementSolution(element, sol, fvGeometry.gridGeometry());

        // resize volume variables to the required size
        volumeVariables_.resize(fvGeometry.numScv());
        for (auto&& scv : scvs(fvGeometry))
            volumeVariables_[scv.indexInElement()].update(elemSol, gridVolVars().problem(), element, scv);
    }

    const VolumeVariables& operator [](std::size_t scvIdx) const
    { return volumeVariables_[scvIdx]; }

    VolumeVariables& operator [](std::size_t scvIdx)
    { return volumeVariables_[scvIdx]; }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return volumeVariables_[scv.indexInElement()]; }

    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    VolumeVariables& operator [](const SubControlVolume& scv)
    { return volumeVariables_[scv.indexInElement()]; }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolVarsPtr_; }

private:
    const GridVolumeVariables* gridVolVarsPtr_;
    std::vector<VolumeVariables> volumeVariables_;
};

} // end namespace Dumux

#endif
