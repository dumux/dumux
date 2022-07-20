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
 * \ingroup CvfeDiscretization
 * \brief The local element solution class for the cvfe method
 */
#ifndef DUMUX_CVFE_ELEMENT_SOLUTION_HH
#define DUMUX_CVFE_ELEMENT_SOLUTION_HH

#include <type_traits>
#include <dune/common/reservedvector.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup CvfeDiscretization
 * \brief The element solution vector
 */
template<class FVElementGeometry, class PV>
class CvfeElementSolution
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    // Dofs are located at corners and in the element center
    static constexpr int numCubeDofs = 1 << dim + 1;

public:
    //! export the primary variables type
    using PrimaryVariables = PV;

    //! Default constructor
    CvfeElementSolution() = default;

    //! Constructor with element and solution and grid geometry
    template<class SolutionVector>
    CvfeElementSolution(const Element& element, const SolutionVector& sol,
                        const GridGeometry& gridGeometry)
    {
        update(element, sol, gridGeometry);
    }

    //! Constructor with element and elemVolVars and fvGeometry
    template<class ElementVolumeVariables>
    CvfeElementSolution(const Element& element, const ElementVolumeVariables& elemVolVars,
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
        const auto numVertices = element.subEntities(GridView::dimension);
        priVars_.resize(numVertices+1);
        for (int vIdx = 0; vIdx < numVertices; ++vIdx)
            priVars_[vIdx] = sol[gridGeometry.dofMapper().subIndex(element, vIdx, GridView::dimension)];

        priVars_[numVertices] = sol[gridGeometry.dofMapper().index(element)];
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

/*!
 * \ingroup CvfeDiscretization
 * \brief  Make an element solution for cvfe schemes
 */
template<class Element, class SolutionVector, class GridGeometry>
auto elementSolution(const Element& element, const SolutionVector& sol, const GridGeometry& gg)
-> std::enable_if_t<GridGeometry::discMethod == DiscretizationMethods::cvfe,
                    CvfeElementSolution<typename GridGeometry::LocalView,
                                      std::decay_t<decltype(std::declval<SolutionVector>()[0])>>
                    >
{
    using PrimaryVariables = std::decay_t<decltype(std::declval<SolutionVector>()[0])>;
    return CvfeElementSolution<typename GridGeometry::LocalView, PrimaryVariables>(element, sol, gg);
}

/*!
 * \ingroup CvfeDiscretization
 * \brief  Make an element solution for cvfe schemes
 */
template<class Element, class ElementVolumeVariables, class FVElementGeometry>
auto elementSolution(const Element& element, const ElementVolumeVariables& elemVolVars, const FVElementGeometry& gg)
-> std::enable_if_t<FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::cvfe,
                    CvfeElementSolution<FVElementGeometry,
                                       typename ElementVolumeVariables::VolumeVariables::PrimaryVariables>>
{
    using PrimaryVariables = typename ElementVolumeVariables::VolumeVariables::PrimaryVariables;
    return CvfeElementSolution<FVElementGeometry, PrimaryVariables>(element, elemVolVars, gg);
}

} // end namespace Dumux

#endif
