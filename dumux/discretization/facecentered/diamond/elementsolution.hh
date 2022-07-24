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

#include <dune/common/reservedvector.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup DiamondDiscretization
 * \brief The global face variables class for face-centered diamond models
 */
template<class FVElementGeometry, class PV>
class FaceCenteredDiamondElementSolution
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int numCubeFaces = 2*dim;

public:
    //! export the primary variables type
    using PrimaryVariables = PV;

    FaceCenteredDiamondElementSolution() = default;

    //! Constructor with element, solution vector and grid geometry
    template<class SolutionVector>
    FaceCenteredDiamondElementSolution(const Element& element, const SolutionVector& sol,
                                       const GridGeometry& gridGeometry)
    {
        update(element, sol, gridGeometry);
    }

    //! Constructor with element, element volume variables and fv element geometry
    template<class ElementVolumeVariables>
    FaceCenteredDiamondElementSolution(const Element& element, const ElementVolumeVariables& elemVolVars,
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
        const auto numFaces = element.subEntities(1);
        priVars_.resize(numFaces);
        for (int fIdx = 0; fIdx < numFaces; ++fIdx)
            priVars_[fIdx] = sol[gridGeometry.dofMapper().subIndex(element, fIdx, 1)];
    }

    //! extract the element solution from the solution vector using a mapper
    template<class SolutionVector>
    void update(const Element& element, const SolutionVector& sol,
                const FVElementGeometry& fvGeometry)
    {
        priVars_.resize(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
            priVars_[scv.localDofIndex()] = sol[scv.dofIndex()];
    }

    //! bracket operator const access
    template<typename IndexType>
    const PrimaryVariables& operator [](IndexType localScvIdx) const
    { return priVars_[localScvIdx]; }

    //! bracket operator
    template<typename IndexType>
    PrimaryVariables& operator [](IndexType localScvIdx)
    { return priVars_[localScvIdx]; }

    //! return the size of the element solution
    std::size_t size() const
    { return priVars_.size(); }

private:
    Dune::ReservedVector<PrimaryVariables, numCubeFaces> priVars_;
};

/*!
 * \ingroup DiamondDiscretization
 * \brief  Make an element solution for face-centered diamond schemes
 */
template<class Element, class SolutionVector, class GridGeometry>
auto elementSolution(const Element& element, const SolutionVector& sol, const GridGeometry& gg)
-> std::enable_if_t<GridGeometry::discMethod == DiscretizationMethods::fcdiamond,
                    FaceCenteredDiamondElementSolution<typename GridGeometry::LocalView,
                                      std::decay_t<decltype(std::declval<SolutionVector>()[0])>>
                    >
{
    using PrimaryVariables = std::decay_t<decltype(std::declval<SolutionVector>()[0])>;
    return FaceCenteredDiamondElementSolution<typename GridGeometry::LocalView, PrimaryVariables>(element, sol, gg);
}

/*!
 * \ingroup DiamondDiscretization
 * \brief  Make an element solution for face-centered diamond schemes
 */
template<class Element, class ElementVolumeVariables, class FVElementGeometry>
auto elementSolution(const Element& element, const ElementVolumeVariables& elemVolVars, const FVElementGeometry& gg)
-> std::enable_if_t<FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::fcdiamond,
                    FaceCenteredDiamondElementSolution<FVElementGeometry, typename ElementVolumeVariables::VolumeVariables::PrimaryVariables>>
{
    using PrimaryVariables = typename ElementVolumeVariables::VolumeVariables::PrimaryVariables;
    return FaceCenteredDiamondElementSolution<FVElementGeometry, PrimaryVariables>(element, elemVolVars, gg);
}

} // end namespace Dumux

#endif
