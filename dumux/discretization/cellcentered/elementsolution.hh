// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCDiscretization
 * \brief The local element solution class for cell-centered methods
 */
#ifndef DUMUX_CC_ELEMENT_SOLUTION_HH
#define DUMUX_CC_ELEMENT_SOLUTION_HH

#include <cassert>
#include <utility>
#include <type_traits>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup CCDiscretization
 * \brief The element solution vector
 */
template<class FVElementGeometry, class PV>
class CCElementSolution
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! export the primary variables type
    using PrimaryVariables = PV;

    //! default constructor
    CCElementSolution() = default;

    //! Constructor with element, solution vector and grid geometry
    template<class SolutionVector>
    CCElementSolution(const Element& element, const SolutionVector& sol,
                      const GridGeometry& gridGeometry)
    : CCElementSolution(sol[gridGeometry.elementMapper().index(element)])
    {}

    //! Constructor with element, element volume variables and fv element geometry
    template<class ElementVolumeVariables>
    CCElementSolution(const Element& element, const ElementVolumeVariables& elemVolVars,
                      const FVElementGeometry& fvGeometry)
    {
        for (const auto& scv : scvs(fvGeometry))
            priVars_ = elemVolVars[scv].priVars();
    }

    //! Constructor with a primary variable object
    CCElementSolution(PrimaryVariables&& priVars)
    : priVars_(std::move(priVars)) {}

    //! Constructor with a primary variable object
    CCElementSolution(const PrimaryVariables& priVars)
    : priVars_(priVars) {}

    //! extract the element solution from the solution vector using a mapper
    template<class SolutionVector>
    void update(const Element& element, const SolutionVector& sol,
                const GridGeometry& gridGeometry)
    {
        priVars_ = sol[gridGeometry.elementMapper().index(element)];
    }

    //! return the size of the element solution
    constexpr std::size_t size() const
    { return 1; }

    //! bracket operator const access
    template<typename IndexType>
    const PrimaryVariables& operator [](IndexType i) const
    {
        assert(i == 0 && "Index exceeds valid range!");
        return priVars_;
    }

    //! bracket operator access
    template<typename IndexType>
    PrimaryVariables& operator [](IndexType i)
    {
        assert(i == 0 && "Index exceeds valid range!");
        return priVars_;
    }

private:
    PrimaryVariables priVars_;
};

/*!
 * \ingroup Discretization
 * \brief  Make an element solution for cell-centered schemes
 */
template<class Element, class SolutionVector, class GridGeometry>
auto elementSolution(const Element& element, const SolutionVector& sol, const GridGeometry& gg)
-> std::enable_if_t<GridGeometry::discMethod == DiscretizationMethods::cctpfa ||
                    GridGeometry::discMethod == DiscretizationMethods::ccmpfa,
                    CCElementSolution<typename GridGeometry::LocalView,
                                      std::decay_t<decltype(std::declval<SolutionVector>()[0])>>
                    >
{
    using PrimaryVariables = std::decay_t<decltype(std::declval<SolutionVector>()[0])>;
    return CCElementSolution<typename GridGeometry::LocalView, PrimaryVariables>(element, sol, gg);
}

/*!
 * \ingroup Discretization
 * \brief  Make an element solution for cell-centered schemes
 */
template<class Element, class ElementVolumeVariables, class FVElementGeometry>
auto elementSolution(const Element& element, const ElementVolumeVariables& elemVolVars, const FVElementGeometry& gg)
-> std::enable_if_t<FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::cctpfa ||
                    FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::ccmpfa,
                    CCElementSolution<FVElementGeometry, typename ElementVolumeVariables::VolumeVariables::PrimaryVariables>>
{
    using PrimaryVariables = typename ElementVolumeVariables::VolumeVariables::PrimaryVariables;
    return CCElementSolution<FVElementGeometry, PrimaryVariables>(element, elemVolVars, gg);
}

/*!
 * \ingroup Discretization
 * \brief  Make an element solution for cell-centered schemes
 * \note This is e.g. used to construct an element solution at Dirichlet boundaries
 */
template<class FVElementGeometry, class PrimaryVariables>
auto elementSolution(PrimaryVariables&& priVars)
-> std::enable_if_t<FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::cctpfa ||
                    FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::ccmpfa,
                    CCElementSolution<FVElementGeometry, PrimaryVariables>>
{
    return CCElementSolution<FVElementGeometry, PrimaryVariables>(std::move(priVars));
}

/*!
 * \ingroup Discretization
 * \brief  Make an element solution for cell-centered schemes
 * \note This is e.g. used to construct an element solution at Dirichlet boundaries
 */
template<class FVElementGeometry, class PrimaryVariables>
auto elementSolution(const PrimaryVariables& priVars)
-> std::enable_if_t<FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::cctpfa ||
                    FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::ccmpfa,
                    CCElementSolution<FVElementGeometry, PrimaryVariables>>
{
    return CCElementSolution<FVElementGeometry, PrimaryVariables>(priVars);
}

} // end namespace Dumux

#endif
