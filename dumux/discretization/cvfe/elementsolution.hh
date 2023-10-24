// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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
#include <dumux/discretization/method.hh>
#include <dumux/discretization/localdoftraits.hh>

namespace Dumux {

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

    // Dofs are located at corners and in the element center
    static constexpr int numCubeDofs = Detail::LocalDofTraits<
        GridView, typename GridGeometry::DiscretizationMethod
    >::numCubeElementDofs;

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
    template<class ElementVolumeVariables>
    CVFEElementSolution(const Element& element, const ElementVolumeVariables& elemVolVars,
                        const FVElementGeometry& fvGeometry)
    {
        priVars_.resize(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
            priVars_[scv.localDofIndex()] = elemVolVars[scv].priVars();
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
                )
            ];
        }
    }

    //! extract the element solution from the solution vector using a local fv geometry
    template<class SolutionVector>
    void update(const Element& element, const SolutionVector& sol,
                const FVElementGeometry& fvGeometry)
    {
        priVars_.resize(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
            priVars_[scv.localDofIndex()] = sol[scv.dofIndex()];
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
    Dune::ReservedVector<PrimaryVariables, numCubeDofs> priVars_;
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
template<class Element, class ElementVolumeVariables, class FVElementGeometry>
auto elementSolution(const Element& element, const ElementVolumeVariables& elemVolVars, const FVElementGeometry& gg)
-> std::enable_if_t<DiscretizationMethods::isCVFE<typename FVElementGeometry::GridGeometry::DiscretizationMethod>,
                    CVFEElementSolution<FVElementGeometry,
                                       typename ElementVolumeVariables::VolumeVariables::PrimaryVariables>>
{
    using PrimaryVariables = typename ElementVolumeVariables::VolumeVariables::PrimaryVariables;
    return CVFEElementSolution<FVElementGeometry, PrimaryVariables>(element, elemVolVars, gg);
}

} // end namespace Dumux

#endif
