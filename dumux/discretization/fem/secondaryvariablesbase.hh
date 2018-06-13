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
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */
#ifndef DUMUX_DISCRETIZATION_FEM_SECONDARY_VARIABLES_BASE_HH
#define DUMUX_DISCRETIZATION_FEM_SECONDARY_VARIABLES_BASE_HH

#include <dumux/implicit/properties.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux
{

/*!
 * \ingroup FemSecondaryVariables
 * \brief Class to store secondary variables for finite element models.
 */
template<class TypeTag>
class FemSecondaryVariablesBase
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    using Element = typename GridView::template Codim<0>::Entity;

public:
    /*!
     * \brief Update all quantities on a given integration point
     *
     * \param elemSol The solution at the dofs connected to this element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param ipData Container holding data on shape functions and gradients at the ip
     */
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const IpData& ipData)
    {
//std::cout << "Wir sind in secVarsBaseUpdate" << std::endl << std::endl;
//printvector(std::cout, elemSol, "secVarsBaseElemSol","");
        // interpolate primary variables
        priVars_ = 0.0;
        for (unsigned int i = 0; i < elemSol.size(); ++i)
        {
        //std::cout << "secVarsBaseElemSolSizeCounter_i: " << i << std::endl;
            PrimaryVariables tmp(elemSol[i]);

            tmp *= ipData.shapeValues()[i];
            priVars_ += tmp;

        //printvector(std::cout, priVars_, "secVarsBasePriVars_","");
        }

        // set the extrusion factor
        extrusionFactor_ = problem.extrusionFactor(element, priVars_);
    }

    /*!
     * \brief Return the vector of primary variables
     */
    const PrimaryVariables& priVars() const
    { return priVars_; }

    /*!
     * \brief Return a component of primary variable vector
     *
     * \param pvIdx The index of the primary variable of interest
     */
    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    /*!
     * \brief Return how much the domain is extruded.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

private:
    PrimaryVariables priVars_;
    Scalar extrusionFactor_;
};

} // end namespace Dumux

#endif
