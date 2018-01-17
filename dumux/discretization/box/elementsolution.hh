// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief The local element solution class for the box method
 */
#ifndef DUMUX_BOX_ELEMENT_SOLUTION_HH
#define DUMUX_BOX_ELEMENT_SOLUTION_HH

#include <dune/istl/bvector.hh>
#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup BoxModel
 * \brief The element solution vector
 */
template<class TypeTag>
class BoxElementSolution
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;

public:
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    //! Default constructors
    BoxElementSolution() = default;

    //! Constructor with element and solution and gridgeometry
    BoxElementSolution(const Element& element, const SolutionVector& sol,
                       const FVGridGeometry& fvGridGeometry)
    {
        update(element, sol, fvGridGeometry);
    }

    //! Constructor with element and solution and elementgeometry (only works for box)
    BoxElementSolution(const Element& element, const SolutionVector& sol,
                       const FVElementGeometry& fvGeometry)
    {
        update(element, sol, fvGeometry);
    }

    //! Constructor with element and elemVolVars and fvGeometry
    BoxElementSolution(const Element& element, const ElementVolumeVariables& elemVolVars,
                       const FVElementGeometry& fvGeometry)
    {
        const auto numVert = element.subEntities(GridView::dimension);
        priVars_.resize(numVert);
        for (const auto& scv : scvs(fvGeometry))
            priVars_[scv.indexInElement()] = elemVolVars[scv].priVars();
    }

    //! extract the element solution from the solution vector using a mapper
    void update(const Element& element, const SolutionVector& sol,
                const FVGridGeometry& fvGridGeometry)
    {
        const auto numVert = element.subEntities(GridView::dimension);
        priVars_.resize(numVert);
        for (int vIdx = 0; vIdx < numVert; ++vIdx)
            priVars_[vIdx] = sol[fvGridGeometry.vertexMapper().subIndex(element, vIdx, GridView::dimension)];
    }

    //! extract the element solution from the solution vector using a local fv geometry
    void update(const Element& element, const SolutionVector& sol,
                const FVElementGeometry& fvGeometry)
    {
        const auto numVert = element.subEntities(GridView::dimension);
        priVars_.resize(numVert);
        for (const auto& scv : scvs(fvGeometry))
            priVars_[scv.indexInElement()] = sol[scv.dofIndex()];
    }

    //! bracket operator const access
    template<typename IndexType>
    const PrimaryVariables& operator [](IndexType i) const
    { return priVars_[i]; }

    //! bracket operator access
    template<typename IndexType>
    PrimaryVariables& operator [](IndexType i)
    { return priVars_[i]; }

private:
    Dune::BlockVector<PrimaryVariables> priVars_;
};

} // end namespace Dumux

#endif
