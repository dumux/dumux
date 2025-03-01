// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief The local variables class
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_ELEMENT_LOCALVARIABLES_HH
#define DUMUX_DISCRETIZATION_CVFE_ELEMENT_LOCALVARIABLES_HH

#include <type_traits>
#include <utility>
#include <vector>

#include <dumux/common/typetraits/localdofs.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/discretization/variablesdeflectionhelper.hh>
#include <dumux/discretization/variablesupdatehelper.hh>

namespace Dumux {

/*!
 * \ingroup CVFEDiscretization
 * \brief The local (stencil) variables class for control-volume finite element
 * \note The class is specialized for versions with and without caching
 * \tparam GLV the grid local variables type
 * \tparam cachingEnabled if the cache is enabled
 */
template<class GLV, bool cachingEnabled>
class CVFEElementLocalVariables;

/*!
 * \ingroup CVFEDiscretization
 * \brief The local (stencil) element variables class for control-volume finite element with caching
 * \note the local variables are stored for the whole grid view in the corresponding GridLocalVariables class
 */
template<class GLV>
class CVFEElementLocalVariables<GLV, /*cachingEnabled*/true>
{
public:
    //! export type of the grid local variables
    using GridLocalVariables = GLV;

    //! export type of the local variables
    using LocalVariables = typename GridLocalVariables::LocalVariables;

    //! export the volume variables type
    using VolumeVariables [[deprecated("Use LocalVariables instead. Will be removed after release 3.10.")]] = LocalVariables;

    //! export type to update local variables
    using UpdateHelper = VariablesUpdateHelper<GridLocalVariables>;

    //! export type to deflect local variables
    template<class FVElementGeometry>
    using DeflectionHelper = VariablesDeflectionHelper<GridLocalVariables, FVElementGeometry>;

    //! Constructor
    CVFEElementLocalVariables(const GridLocalVariables& gridLocalVars)
    : gridLocalVarsPtr_(&gridLocalVars) {}

    const LocalVariables& operator [](std::size_t localDofIdx) const
    { return gridLocalVars().localVars(eIdx_, localDofIdx); }

    template<class ScvOrLocalDof, typename std::enable_if_t<!std::is_integral<ScvOrLocalDof>::value, int> = 0>
    const LocalVariables& operator [](const ScvOrLocalDof& scvOrLocalDof) const
    {
        if constexpr (Detail::isLocalDofType<ScvOrLocalDof>())
            return gridLocalVars().localVars(eIdx_, scvOrLocalDof.index());
        else
            return gridLocalVars().localVars(eIdx_, scvOrLocalDof.indexInElement());
    }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    CVFEElementLocalVariables bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
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
    CVFEElementLocalVariables bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
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

    //! The grid local variables object we are a restriction of
    const GridLocalVariables& gridLocalVars() const
    { return *gridLocalVarsPtr_; }

private:
    const GridLocalVariables* gridLocalVarsPtr_;
    std::size_t eIdx_;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief The local (stencil) element variables class for control-volume finite element without caching
 */
template<class GLV>
class CVFEElementLocalVariables<GLV, /*cachingEnabled*/false>
{
public:
    //! export type of the grid local variables
    using GridLocalVariables = GLV;

    //! export type of the local variables
    using LocalVariables = typename GridLocalVariables::LocalVariables;

    //! export the volume variables type
    using VolumeVariables [[deprecated("Use LocalVariables instead. Will be removed after release 3.10.")]] = LocalVariables;

    //! export type to update local variables
    using UpdateHelper = VariablesUpdateHelper<GridLocalVariables>;

    //! export type to deflect local variables
    template<class FVElementGeometry>
    using DeflectionHelper = VariablesDeflectionHelper<GridLocalVariables, FVElementGeometry>;

    //! Constructor
    CVFEElementLocalVariables(const GridLocalVariables& gridLocalVars)
    : gridLocalVarsPtr_(&gridLocalVars) {}

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    template<class FVElementGeometry, class SolutionVector>
    CVFEElementLocalVariables bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
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
    CVFEElementLocalVariables bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
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
        localVariables_.resize(Detail::numLocalDofs(fvGeometry));

        // update cv related dofs where there exists a localDof
        for (const auto& localDof : localDofs(fvGeometry))
            localVariables_[localDof.index()].update(elemSol, gridLocalVars().problem(), fvGeometry, localDof);
    }

    const LocalVariables& operator [](std::size_t localIdx) const
    { return localVariables_[localIdx]; }

    LocalVariables& operator [](std::size_t localIdx)
    { return localVariables_[localIdx]; }

    template<class ScvOrLocalDof, typename std::enable_if_t<!std::is_integral<ScvOrLocalDof>::value, int> = 0>
    const LocalVariables& operator [](const ScvOrLocalDof& scvOrLocalDof) const
    {
        if constexpr (Detail::isLocalDofType<ScvOrLocalDof>())
            return localVariables_[scvOrLocalDof.index()];
        else
            return localVariables_[scvOrLocalDof.indexInElement()];
    }

    template<class ScvOrLocalDof, typename std::enable_if_t<!std::is_integral<ScvOrLocalDof>::value, int> = 0>
    LocalVariables& operator [](const ScvOrLocalDof& scvOrLocalDof)
    {
        if constexpr (Detail::isLocalDofType<ScvOrLocalDof>())
            return localVariables_[scvOrLocalDof.index()];
        else
            return localVariables_[scvOrLocalDof.indexInElement()];
    }

    //! The grid local variables object we are a restriction of
    const GridLocalVariables& gridLocalVars() const
    { return *gridLocalVarsPtr_; }

private:
    const GridLocalVariables* gridLocalVarsPtr_;
    std::vector<LocalVariables> localVariables_;
};

} // end namespace Dumux

#endif
