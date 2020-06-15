// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredElementSolution
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_ELEMENT_SOLUTION_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_ELEMENT_SOLUTION_HH



#include <type_traits>
#include <array>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief The global face variables class for staggered models
 */
template<class FVElementGeometry, class PV>
class FaceCenteredStaggeredElementSolution
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;

    static constexpr auto numScvsPerElement = GridView::Grid::dimension * 2;

public:
    //! export the primary variables type
    using PrimaryVariables = PV;

    FaceCenteredStaggeredElementSolution() = default;

    //! Constructor with element, solution vector and grid geometry
    template<class Element, class SolutionVector>
    FaceCenteredStaggeredElementSolution(const Element& element,
                                         const SolutionVector& sol,
                                         const GridGeometry& gridGeometry)
    {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(element);

        for (const auto& scv : scvs(fvGeometry))
            priVars_[scv.indexInElement()] = sol[scv.dofIndex()];
    }

    //! Constructor with element, element volume variables and fv element geometry
    template<class ElementVolumeVariables>
    FaceCenteredStaggeredElementSolution(const Element& element,
                                         const ElementVolumeVariables& elemVolVars,
                                         const FVElementGeometry& fvGeometry)
    {
        for (const auto& scv : scvs(fvGeometry))
            priVars_ = elemVolVars[scv].priVars();
    }

    //! Constructor with a primary variable object
    FaceCenteredStaggeredElementSolution(PrimaryVariables&& priVars)
    {
        priVars_[0] = std::move(priVars);
        for (int i = 1; i < priVars_.size(); ++i)
           priVars_[i] = priVars_[0];
    }

    //! Constructor with a primary variable object
    FaceCenteredStaggeredElementSolution(const PrimaryVariables& priVars)
    {
        priVars_[0] = priVars;
        for (int i = 1; i < priVars_.size(); ++i)
           priVars_[i] = priVars_[0];
    }

    //! extract the element solution from the solution vector using a mapper
    template<class SolutionVector>
    void update(const Element& element, const SolutionVector& sol,
                const GridGeometry& gridGeometry)
    {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(element);

        for (const auto& scv : scvs(fvGeometry))
            priVars_[scv.indexInElement()] = sol[scv.dofIndex()];
    }

    //! bracket operator const access
    const PrimaryVariables& operator [](SmallLocalIndexType localScvIdx) const
    {
        return priVars_[localScvIdx];
    }

    //! bracket operator
    PrimaryVariables& operator [](SmallLocalIndexType localScvIdx)
    {
         return priVars_[localScvIdx];
    }

    //! return the size of the element solution
    constexpr std::size_t size() const
    { return numScvsPerElement; }

private:

    std::array<PrimaryVariables, numScvsPerElement> priVars_;
};

/*!
 * \ingroup CCDiscretization
 * \brief  Make an element solution for face-centered staggered schemes
 */
template<class Element, class SolutionVector, class GridGeometry>
auto elementSolution(const Element& element, const SolutionVector& sol, const GridGeometry& gg)
-> std::enable_if_t<GridGeometry::discMethod == DiscretizationMethod::fcstaggered,
                    FaceCenteredStaggeredElementSolution<typename GridGeometry::LocalView,
                                      std::decay_t<decltype(std::declval<SolutionVector>()[0])>>
                    >
{
    using PrimaryVariables = std::decay_t<decltype(std::declval<SolutionVector>()[0])>;
    return FaceCenteredStaggeredElementSolution<typename GridGeometry::LocalView, PrimaryVariables>(element, sol, gg);
}

/*!
 * \ingroup CCDiscretization
 * \brief  Make an element solution for face-centered staggered schemes
 */
template<class Element, class ElementVolumeVariables, class FVElementGeometry>
auto elementSolution(const Element& element, const ElementVolumeVariables& elemVolVars, const FVElementGeometry& gg)
-> std::enable_if_t<FVElementGeometry::GridGeometry::discMethod == DiscretizationMethod::fcstaggered,
                    FaceCenteredStaggeredElementSolution<FVElementGeometry, typename ElementVolumeVariables::VolumeVariables::PrimaryVariables>>
{
    using PrimaryVariables = typename ElementVolumeVariables::VolumeVariables::PrimaryVariables;
    return FaceCenteredStaggeredElementSolution<FVElementGeometry, PrimaryVariables>(element, elemVolVars, gg);
}

/*!
 * \ingroup CCDiscretization
 * \brief  Make an element solution for face-centered staggered schemes
 * \note This is e.g. used to contruct an element solution at Dirichlet boundaries
 */
template<class FVElementGeometry, class PrimaryVariables>
auto elementSolution(PrimaryVariables&& priVars)
-> std::enable_if_t<FVElementGeometry::GridGeometry::discMethod == DiscretizationMethod::fcstaggered,
                    FaceCenteredStaggeredElementSolution<FVElementGeometry, std::decay_t<PrimaryVariables>>>
{
    return FaceCenteredStaggeredElementSolution<FVElementGeometry, std::decay_t<PrimaryVariables>>(std::forward<PrimaryVariables>(priVars));
}

} // end namespace Dumux

#endif
