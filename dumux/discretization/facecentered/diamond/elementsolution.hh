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
 * \ingroup DiamondDiscretization
 * \copydoc Dumux::FaceCenteredDiamondElementSolution
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_ELEMENT_SOLUTION_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_ELEMENT_SOLUTION_HH



#include <type_traits>
#include <array>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup DiamondDiscretization
 * \brief The global face variables class for staggered models
 */
template<class FVElementGeometry, class PV>
class FaceCenteredDiamondElementSolution
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using LocalIndexType = typename FVElementGeometry::SubControlVolume::LocalIndexType;

public:
    //! export the primary variables type
    using PrimaryVariables = PV;

    FaceCenteredDiamondElementSolution() = default;

    //! Constructor with element, solution vector and grid geometry
    template<class Element, class SolutionVector>
    FaceCenteredDiamondElementSolution(const Element& element,
                                       const SolutionVector& sol,
                                       const GridGeometry& gridGeometry)
    {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(element);

        priVars_.resize(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
            priVars_[scv.indexInElement()] = sol[scv.dofIndex()];
    }

    //! Constructor with element, element volume variables and fv element geometry
    template<class ElementVolumeVariables>
    FaceCenteredDiamondElementSolution(const Element& element,
                                       const ElementVolumeVariables& elemVolVars,
                                       const FVElementGeometry& fvGeometry)
    {
        priVars_.resize(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
            priVars_[scv.indexInElement()] = elemVolVars[scv].priVars();
    }

    //! Constructor with a primary variable object
    FaceCenteredDiamondElementSolution(const Element& element,
                                       PrimaryVariables&& priVars,
                                       const FVElementGeometry& fvGeometry)
    {
        priVars_.resize(fvGeometry.numScv());
        priVars_[0] = std::move(priVars);
        for (int i = 1; i < priVars_.size(); ++i)
           priVars_[i] = priVars_[0];
    }

    //! Constructor with a primary variable object
    FaceCenteredDiamondElementSolution(const Element& element,
                                       const PrimaryVariables& priVars,
                                       const FVElementGeometry& fvGeometry)
    {
        priVars_.resize(fvGeometry.numScv());
        priVars_[0] = priVars;
        for (int i = 1; i < priVars_.size(); ++i)
           priVars_[i] = priVars_[0];
    }

    //! extract the element solution from the solution vector using a mapper
    template<class SolutionVector>
    void update(const Element& element,
                const SolutionVector& sol,
                const GridGeometry& gridGeometry)
    {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(element);

        priVars_.resize(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
            priVars_[scv.indexInElement()] = sol[scv.dofIndex()];
    }

    //! bracket operator const access
    const PrimaryVariables& operator [](LocalIndexType localScvIdx) const
    {
        return priVars_[localScvIdx];
    }

    //! bracket operator
    PrimaryVariables& operator [](LocalIndexType localScvIdx)
    {
         return priVars_[localScvIdx];
    }

    //! return the size of the element solution
    std::size_t size() const
    { return priVars_.size(); }

private:
    std::vector<PrimaryVariables> priVars_;
};

/*!
 * \ingroup CCDiscretization
 * \brief  Make an element solution for face-centered staggered schemes
 */
template<class Element, class SolutionVector, class GridGeometry>
auto elementSolution(const Element& element, const SolutionVector& sol, const GridGeometry& gg)
-> std::enable_if_t<GridGeometry::discMethod == DiscretizationMethod::fcdiamond,
                    FaceCenteredDiamondElementSolution<typename GridGeometry::LocalView,
                                      std::decay_t<decltype(std::declval<SolutionVector>()[0])>>
                    >
{
    using PrimaryVariables = std::decay_t<decltype(std::declval<SolutionVector>()[0])>;
    return FaceCenteredDiamondElementSolution<typename GridGeometry::LocalView, PrimaryVariables>(element, sol, gg);
}

/*!
 * \ingroup CCDiscretization
 * \brief  Make an element solution for face-centered staggered schemes
 */
template<class Element, class ElementVolumeVariables, class FVElementGeometry>
auto elementSolution(const Element& element, const ElementVolumeVariables& elemVolVars, const FVElementGeometry& gg)
-> std::enable_if_t<FVElementGeometry::GridGeometry::discMethod == DiscretizationMethod::fcdiamond,
                    FaceCenteredDiamondElementSolution<FVElementGeometry, typename ElementVolumeVariables::VolumeVariables::PrimaryVariables>>
{
    using PrimaryVariables = typename ElementVolumeVariables::VolumeVariables::PrimaryVariables;
    return FaceCenteredDiamondElementSolution<FVElementGeometry, PrimaryVariables>(element, elemVolVars, gg);
}

/*!
 * \ingroup CCDiscretization
 * \brief  Make an element solution for face-centered staggered schemes
 * \note This is e.g. used to contruct an element solution at Dirichlet boundaries
 */
template<class FVElementGeometry, class PrimaryVariables>
auto elementSolution(PrimaryVariables&& priVars)
-> std::enable_if_t<FVElementGeometry::GridGeometry::discMethod == DiscretizationMethod::fcdiamond,
                    FaceCenteredDiamondElementSolution<FVElementGeometry, std::decay_t<PrimaryVariables>>>
{
    return FaceCenteredDiamondElementSolution<FVElementGeometry, std::decay_t<PrimaryVariables>>(std::forward<PrimaryVariables>(priVars));
}

} // end namespace Dumux

#endif
