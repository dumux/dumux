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
 * \brief Quantities required by the linear elasticity box
 *        model defined on a vertex.
 */


#ifndef DUMUX_ELASTIC_SECONDARY_VARIABLES_HH
#define DUMUX_ELASTIC_SECONDARY_VARIABLES_HH

#include <dumux/discretization/fem/secondaryvariablesbase.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup ElasticFemModel
 * \ingroup FemImplicitSecondaryVariables
 * \brief Contains the quantities which are constant within a
 *        finite volume in the linear elasticity model.
 */
template <class TypeTag>
class ElasticSecondaryVariables : public FemSecondaryVariablesBase<TypeTag>
{
    using ParentType = FemSecondaryVariablesBase<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MechanicalLaw = typename GET_PROP_TYPE(TypeTag, MechanicalLaw);
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static constexpr int dim = GridView::dimension;
    using DimVector = Dune::FieldVector<Scalar, dim>;

    using Element = typename GridView::template Codim<0>::Entity;

public:
    /*!
     * \copydoc SecondaryVariablesBase::update
     */
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const IpData& ipData)
    { //std::cout << "KAMEN HIER VORBEI secvars" << std::endl;
        ParentType::update(elemSol, problem, element, ipData);

        for (int i = 0; i < dim; ++i)
            displacement_[i] = this->priVar(Indices::u(i));

        // retrieve Lame parameters and rock density from spatialParams
        const auto& lameParams = problem.spatialParams().lameParams(element, this->priVars());
        lambda_ = lameParams[0];
        mu_ = lameParams[1];

        // the density of the solid material
        rockDensity_ = problem.spatialParams().rockDensity(element, this->priVars());
    }

    /*!
      * \brief Return the Lame parameter lambda \f$\mathrm{[Pa]}\f$ at the integration point.
      */
    Scalar lambda() const
    { return lambda_; }

    /*!
      * \brief Return the Lame parameter mu \f$\mathrm{[Pa]}\f$ at the integration point.
      */
    Scalar mu() const
    { return mu_; }

    /*!
     * \brief Returns the rock density \f$\mathrm{[kg / m^3]}\f$ at the integration point.
     */
    Scalar rockDensity() const
    { return rockDensity_; }

    /*!
     * \brief Returns the solid displacement \f$\mathrm{[m]}\f$ in space
     * directions dimIdx at the integration point.
     */
    Scalar displacement(int dimIdx) const
    { return displacement_[dimIdx]; }

    /*!
     * \brief Returns the solid displacement vector \f$\mathrm{[m]}\f$
     *  at the integration point.
     */
    const DimVector& displacement() const
    { return displacement_; }

protected:
    DimVector displacement_;
    Scalar lambda_;
    Scalar mu_;
    Scalar rockDensity_;
};

}

#endif
