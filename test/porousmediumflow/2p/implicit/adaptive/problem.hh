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
 * \ingroup TwoPTests
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_LENSPROBLEM_ADAPTIVE_HH
#define DUMUX_LENSPROBLEM_ADAPTIVE_HH

#include "../incompressible/problem.hh"

#include <dumux/io/container.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */
template <class TypeTag >
class TwoPTestProblemAdaptive : public TwoPTestProblem<TypeTag>
{
    using ParentType = TwoPTestProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimensionworld>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethod::box;

public:
    TwoPTestProblemAdaptive(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        if(!isBox)
            initialValues_ = readFileToContainer<std::vector<PrimaryVariables>>("initialsolutioncc.txt");
        else
            initialValues_ = readFileToContainer<std::vector<PrimaryVariables>>("initialsolutionbox.txt");
    }

    /*!
     * \brief Evaluates the initial value for an element for cell-centered models.
     */
    PrimaryVariables initial(const Element& element) const
    {
        const auto delta = 0.0625;
        unsigned int cellsX = this->gridGeometry().bBoxMax()[0]/delta;
        const auto globalPos = element.geometry().center();

        // the input data corresponds to a uniform grid with discretization length deltaX_
        unsigned int dataIdx = std::trunc(globalPos[1]/delta) * cellsX + std::trunc(globalPos[0]/delta);
        return initialValues_[dataIdx];
    }

    /*!
     * \brief Evaluates the initial value for a vertex for vertex-centered models.
     */
    PrimaryVariables initial(const Vertex& vertex) const
    {
        const auto delta = 0.0625;
        unsigned int verticesX = this->gridGeometry().bBoxMax()[0]/delta + 1;
        const auto globalPos = vertex.geometry().center();

        // the input data corresponds to a uniform grid with discretization length deltaX_
        unsigned int dataIdx = std::trunc(globalPos[1]/delta) * verticesX + std::trunc(globalPos[0]/delta);
        return initialValues_[dataIdx];
    }

private:
    std::vector<PrimaryVariables> initialValues_;
};

} // end namespace Dumux

#endif
