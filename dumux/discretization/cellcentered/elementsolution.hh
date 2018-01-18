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
 * \ingroup CCDiscretization
 * \brief The local element solution class for cell-centered methods
 */
#ifndef DUMUX_CC_ELEMENT_SOLUTION_HH
#define DUMUX_CC_ELEMENT_SOLUTION_HH

#include <dune/istl/bvector.hh>
#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup CCDiscretization
 * \brief The element solution vector
 */
template<class TypeTag>
class CCElementSolution
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

public:
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    CCElementSolution() = default;

    //! Constructor with element and solution
    CCElementSolution(const Element& element, const SolutionVector& sol,
                      const FVGridGeometry& fvGridGeometry)
    : CCElementSolution(sol[fvGridGeometry.elementMapper().index(element)])
    {}

    //! Constructor with element and elemVolVars and fvGeometry
    CCElementSolution(const Element& element, const ElementVolumeVariables& elemVolVars,
                      const FVElementGeometry& fvGeometry)
    {
        for(const auto& scv : scvs(fvGeometry))
            priVars_ = elemVolVars[scv].priVars();
    }

    //! Constructor with a primary variable object
    CCElementSolution(PrimaryVariables&& priVars)
    : priVars_(std::move(priVars)) {}

    //! Constructor with a primary variable object
    CCElementSolution(const PrimaryVariables& priVars)
    : priVars_(priVars) {}

    //! extract the element solution from the solution vector using a mapper
    void update(const Element& element, const SolutionVector& sol,
                const FVGridGeometry& fvGridGeometry)
    {
        priVars_ = sol[fvGridGeometry.elementMapper().index(element)];
    }

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

} // end namespace Dumux

#endif
