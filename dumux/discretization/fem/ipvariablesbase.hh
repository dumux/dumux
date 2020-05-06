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
 * \ingroup FEMDiscretization
 * \brief Base class for the integration point variables that store
 *        the primary and secondary variables evaluated at an integration point.
 *        This base class solely stores the primary variables, derive from this
 *        class to implement primary variables-dependent secondary variables being
 *        evalated at integration points.
 */
#ifndef DUMUX_FEM_INTEGRATIONPOINT_VARIABLES_BASE_HH
#define DUMUX_FEM_INTEGRATIONPOINT_VARIABLES_BASE_HH

namespace Dumux {

/*!
 * \ingroup FEMDiscretization
 * \brief Base class for the integration point variables that store
 *        the primary and secondary variables evaluated at an integration point.
 *        This base class solely stores the primary variables, derive from this
 *        class to implement primary variables-dependent secondary variables being
 *        evalated at integration points.
 * \tparam PV The type used for primary variables
 */
template<class PV>
class IntegrationPointVariablesBase
{
    using Scalar = typename PV::value_type;

public:
    //! export the type used for the primary variables
    using PrimaryVariables = PV;

    /*!
     * \brief Update all quantities for a given integration point
     * \param elemSol The solution at the dofs connected to this element
     * \param problem The problem to be solved
     * \param element The element
     * \param ipData Container holding data on shape functions and gradients
     *               at the integration point
     */
    template<class ElemSol, class Problem, class Element, class IpData>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const IpData& ipData)
    {
        // interpolate primary variables
        priVars_ = 0.0;
        for (unsigned int i = 0; i < elemSol.size(); ++i)
        {
            PrimaryVariables tmp(elemSol[i]);
            tmp *= ipData.shapeValues()[i];
            priVars_ += tmp;
        }

        extrusionFactor_ = problem.extrusionFactor(element, ipData, elemSol);
    }

    //! Return a component of primary variable vector for a given index
    Scalar priVar(unsigned int pvIdx) const
    { return priVars_[pvIdx]; }

    //! Return the vector of primary variables
    const PrimaryVariables& priVars() const
    { return priVars_; }

    //! Extrusion of the domain
    const Scalar extrusionFactor() const
    { return extrusionFactor_; }

private:
    // data members
    Scalar extrusionFactor_;
    PrimaryVariables priVars_;
};

} // end namespace Dumux

#endif
