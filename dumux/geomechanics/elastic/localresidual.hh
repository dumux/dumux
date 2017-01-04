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
 *
 * \brief Element-wise calculation the local Jacobian for the
 *        linear elastic model in the fully implicit scheme.
 */
#ifndef DUMUX_ELASTIC_LOCAL_RESIDUAL_HH
#define DUMUX_ELASTIC_LOCAL_RESIDUAL_HH

#include "properties.hh"

namespace Dumux
{
/*!
 *
 * \ingroup ElasticFemModel
 * \ingroup FemImplicitLocalResidual
 * \brief Calculate the local Jacobian for the linear
 *        elasticity model
 *
 * This class is used to fill the gaps in FemLocalResidual for
 * the linear elasticity model.
 */
template<class TypeTag>
class ElasticLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseLocalResidual);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MechanicalLaw = typename GET_PROP_TYPE(TypeTag, MechanicalLaw);
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SecondaryVariables = typename GET_PROP_TYPE(TypeTag, SecondaryVariables);

    static constexpr int dimWorld = GridView::dimension;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    using typename ParentType::FluxTermType;

    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        within a finite volume.
     *
     *        \param element The finite element
     *        \param ipData Data on shape values and gradients at the integration point
     *        \param secVars Secondary variables object evaluated at integration point
     *        \param elemSol The current primary variables at the dofs of the element
     */
    PrimaryVariables computeStorage(const Element& element,
                                    const IpData& ipData,
                                    const SecondaryVariables& secVars,
                                    const ElementSolution& elemSol) const
    {
        // quasistationary conditions assumed
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Evaluate the stresses.
     *
     *        \param element The finite element
     *        \param ipData Data on shape values and gradients at the integration point
     *        \param secVars Secondary variables object evaluated at integration point
     *        \param elemSol The current primary variables at the dofs of the element
     */
    FluxTermType computeFlux(const Element& element,
                             const IpData& ipData,
                             const SecondaryVariables& secVars,
                             const ElementSolution& elemSol) const
    {
        const auto& lameParams = this->problem().spatialParams().lameParams(element, secVars.priVars());
        return MechanicalLaw::stressTensor(element, ipData, secVars, elemSol, lameParams);
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     *        \param element The finite element
     *        \param ipData Data on shape values and gradients at the integration point
     *        \param secVars Secondary variables object evaluated at integration point
     *        \param elemSol The current primary variables at the dofs of the element
     *
     */
    PrimaryVariables computeSource(const Element& element,
                                   const IpData& ipData,
                                   const SecondaryVariables& secVars,
                                   const ElementSolution& elemSol) const
    {
        PrimaryVariables source(0.0);

        source += ParentType::computeSource(element, ipData, secVars);

        // gravity term of the solid matrix in the momentum balance
        GlobalPosition gravityTerm(0.0);
        gravityTerm = this->problem().gravityAtPos(ipData.ipGlobal());
        gravityTerm *= secVars.rockDensity();

        for (int i = 0; i < dimWorld; ++i)
          source[Indices::momentum(i)] += gravityTerm[i];

        return source;
    }

};

} // end namespace Dumux

#endif // DUMUX_ELASTIC_LOCAL_RESIDUAL_HH
