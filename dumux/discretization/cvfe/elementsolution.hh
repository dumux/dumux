// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief The local element solution class for control-volume finite element methods
 */
#ifndef DUMUX_CVFE_ELEMENT_SOLUTION_HH
#define DUMUX_CVFE_ELEMENT_SOLUTION_HH

#include <type_traits>
#include <dune/common/reservedvector.hh>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux {

namespace Detail {

template<class ElementVariables>
consteval auto primaryVariablesType()
{
    if constexpr (requires { typename ElementVariables::Variables::PrimaryVariables; })
        return std::type_identity<typename ElementVariables::Variables::PrimaryVariables>{};
    else if constexpr (requires { typename ElementVariables::VolumeVariables::PrimaryVariables; })
        return std::type_identity<typename ElementVariables::VolumeVariables::PrimaryVariables>{};
    else
        return std::type_identity<void>{};
}

template<class ElementVariables>
using PrimaryVariables_t = typename decltype(primaryVariablesType<ElementVariables>())::type;

} // end namespace Detail

/*!
 * \ingroup CVFEDiscretization
 * \brief The element solution vector
 */
template<class FVElementGeometry, class PV>
class CVFEElementSolution
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;

    // Maximum number of local element dofs
    static constexpr int maxNumLocalDofs = Detail::LocalDofs::maxNumLocalDofs<FVElementGeometry>();

public:
    //! export the primary variables type
    using PrimaryVariables = PV;

    //! Default constructor
    CVFEElementSolution() = default;

    //! Constructor with element and solution and grid geometry
    template<class SolutionVector>
    CVFEElementSolution(const Element& element, const SolutionVector& sol,
                        const GridGeometry& gridGeometry)
    {
        update(element, sol, gridGeometry);
    }

    //! Constructor with element and elemVolVars and fvGeometry
    template<class ElementVariables>
    CVFEElementSolution(const Element& element, const ElementVariables& elemVars,
                        const FVElementGeometry& fvGeometry)
    {
        priVars_.resize(Detail::LocalDofs::numLocalDofs(fvGeometry));
        for (const auto& localDof : localDofs(fvGeometry))
            priVars_[localDof.index()] = elemVars[localDof.index()].priVars();
    }

    //! extract the element solution from the solution vector using a mapper
    template<class SolutionVector>
    void update(const Element& element, const SolutionVector& sol,
                const GridGeometry& gridGeometry)
    {
        // TODO: this implementation only works if there is only one dof per codim/entity
        // As local-global mappings are provided by the grid geometry
        // this logic should maybe become part of the grid geometry.
        const auto& localCoeff = gridGeometry.feCache().get(element.type()).localCoefficients();
        priVars_.resize(localCoeff.size());
        for (int localDofIdx = 0; localDofIdx < localCoeff.size(); ++localDofIdx)
        {
            const auto& localKey = localCoeff.localKey(localDofIdx);
            priVars_[localDofIdx] = sol[
                gridGeometry.dofMapper().subIndex(
                    element, localKey.subEntity(), localKey.codim()
                ) + localKey.index()
            ];
        }
    }

    //! extract the element solution from the solution vector using a local fv geometry
    template<class SolutionVector>
    void update(const Element& element, const SolutionVector& sol,
                const FVElementGeometry& fvGeometry)
    {
        priVars_.resize(Detail::LocalDofs::numLocalDofs(fvGeometry));
        for (const auto& localDof : localDofs(fvGeometry))
            priVars_[localDof.index()] = sol[localDof.dofIndex()];
    }

    //! return the size of the element solution
    std::size_t size() const
    { return priVars_.size(); }

    //! bracket operator const access
    template<typename IndexType>
    const PrimaryVariables& operator [](IndexType i) const
    { return priVars_[i]; }

    //! bracket operator access
    template<typename IndexType>
    PrimaryVariables& operator [](IndexType i)
    { return priVars_[i]; }

private:
    Dune::ReservedVector<PrimaryVariables, maxNumLocalDofs> priVars_;
};

//! Make an element solution for control-volume finite element schemes
template<class Element, class SolutionVector, class GridGeometry>
auto elementSolution(const Element& element, const SolutionVector& sol, const GridGeometry& gg)
-> std::enable_if_t<DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>,
                    CVFEElementSolution<typename GridGeometry::LocalView,
                                      std::decay_t<decltype(std::declval<SolutionVector>()[0])>>
                    >
{
    using PrimaryVariables = std::decay_t<decltype(std::declval<SolutionVector>()[0])>;
    return CVFEElementSolution<typename GridGeometry::LocalView, PrimaryVariables>(element, sol, gg);
}

//!  Make an element solution for control-volume finite element schemes
template<class Element, class ElementVariables, class FVElementGeometry>
auto elementSolution(const Element& element, const ElementVariables& elemVars, const FVElementGeometry& gg)
-> std::enable_if_t<DiscretizationMethods::isCVFE<typename FVElementGeometry::GridGeometry::DiscretizationMethod>
                    && !std::is_void_v<Detail::PrimaryVariables_t<ElementVariables>>,
                    CVFEElementSolution<FVElementGeometry,
                                        Detail::PrimaryVariables_t<ElementVariables>>
                    >
{
    using PrimaryVariables = Detail::PrimaryVariables_t<ElementVariables>;
    return CVFEElementSolution<FVElementGeometry, PrimaryVariables>(element, elemVars, gg);
}

} // end namespace Dumux

#endif
