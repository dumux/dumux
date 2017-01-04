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
 * \brief Base class for all models which use the linear elasticity model.
 *        Adaption of the fully implicit scheme to the linear elasticity model.
 */
#ifndef DUMUX_ELASTIC_MODEL_HH
#define DUMUX_ELASTIC_MODEL_HH

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup ElasticBoxModel
 * \brief Adaption of the fully implicit scheme to the linear elasticity model.
 *
 * This model implements a linear elastic solid using Hooke's law as
 * stress-strain relation and a quasi-stationary momentum balance equation:
 \f[
  \boldsymbol{\sigma} = 2\,G\,\boldsymbol{\epsilon} + \lambda \,\text{tr} (\boldsymbol{\epsilon}) \, \boldsymbol{I}.
 \f]
 *
 * with the strain tensor \f$\boldsymbol{\epsilon}\f$ as a function of the solid displacement gradient \f$\textbf{grad} \boldsymbol{u}\f$:
 \f[
  \boldsymbol{\epsilon} = \frac{1}{2} \, (\textbf{grad} \boldsymbol{u} + \textbf{grad}^T \boldsymbol{u}).
 \f]
 *
 * Gravity can be enabled or disabled via the property system.
 * By inserting this into the momentum balance equation, one gets
 \f[
 \text{div} \boldsymbol{\sigma} + \varrho {\textbf g} = 0 \;,
 \f]
 *
 * The equation is discretized using a vertex-centered finite volume (box)
 * scheme as spatial discretization.
 *
 */

template<class TypeTag >
class ElasticModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MechanicalLaw = typename GET_PROP_TYPE(TypeTag, MechanicalLaw);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SecondaryVariables = typename GET_PROP_TYPE(TypeTag, SecondaryVariables);

    static const int dim = GridView::dimension;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

public:
    /*!
     * \brief \copybrief ImplicitModel::addOutputVtkFields
     *
     * Specialization for the ElasticBoxModel, adding solid displacement,
     * stresses and the process rank to the VTK writer.
     */
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol, MultiWriter &writer)
    {
        // create the required scalar fields
        unsigned numVert = this->gridView_().size(dim);
        unsigned numElements = this->gridView_().size(0);

        auto& ux = *writer.allocateManagedBuffer(numVert);
        auto& uy = *writer.allocateManagedBuffer(numVert);
        auto& uz = *writer.allocateManagedBuffer(numVert);
        auto& sigmax = *writer.template allocateManagedBuffer<Scalar, dim>(numElements);
        auto& sigmay = *writer.template allocateManagedBuffer<Scalar, dim>(numElements);
        auto& sigmaz = *writer.template allocateManagedBuffer<Scalar, dim>(numElements);
        auto& rank = *writer.allocateManagedBuffer(numElements);

        // initialize stress fields
        for (unsigned int i = 0; i < numElements; ++i)
        {
            sigmax[i] = 0;
            if (dim > 1)
                sigmay[i] = 0;
            if (dim > 2)
                sigmaz[i] = 0;
        }

        auto localView = this->feBasis().localView();
        auto localIndexSet = this->feBasis().localIndexSet();

        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
        {
            auto eIdx = this->problem_().model().elementMapper().index(element);

            // rank output
            rank[eIdx] = this->gridView_().comm().rank();

            // bind local restrictions to element
            localView.bind(element);
            localIndexSet.bind(localView);

            // obtain the finite element
            const auto& fe = localView.tree().finiteElement();

            // loop over dofs inside the element and store solution at vertices
            // TODO HOW TO INCLUDE SUBSAMPLING
            const auto numLocalDofs = fe.localBasis().size();
            ElementSolution elemSol(numLocalDofs);
            for (int i = 0; i < numLocalDofs; ++i)
            {
                // only proceed for vertex dofs
                if (fe.localCoefficients().localKey(i).codim() != dim)
                    continue;

                auto dofIdxGlobal = localIndexSet.index(i);

                auto dofSol = sol[dofIdxGlobal];
                ux[dofIdxGlobal] = dofSol[Indices::u(0)];
                if (dim >= 2)
                    uy[dofIdxGlobal] = dofSol[Indices::u(1)];
                if (dim >= 3)
                    uz[dofIdxGlobal] = dofSol[Indices::u(2)];

                elemSol[i] = std::move(dofSol);
            }

            // obtain element geometry
            auto eg = element.geometry();

            // evaluate shape function data and secondary variables at the cell center
            IpData ipData(eg, eg.local(eg.center()), fe.localBasis());
            SecondaryVariables secVars;
            secVars.update(elemSol, this->problem_(), element, ipData);

            // get the lame parameters
            const auto& lameParams = this->problem_().spatialParams().lameParams(element, secVars.priVars());

            // compute the stress tensor and add to container
            auto sigma = MechanicalLaw::stressTensor(element, ipData, secVars, elemSol, lameParams);
            sigmax[eIdx] += sigma[0];
            if (dim >= 2)
                sigmay[eIdx] += sigma[1];
            if (dim == 3)
                sigmaz[eIdx] += sigma[2];
        }

        writer.attachDofData(ux, "ux", true);
        if (dim >= 2)
            writer.attachDofData(uy, "uy", true);
        if (dim == 3)
            writer.attachDofData(uz, "uz", true);
        writer.attachCellData(sigmax, "stress X", dim);
        if (dim >= 2)
        writer.attachCellData(sigmay, "stress Y", dim);
        if (dim == 3)
        writer.attachCellData(sigmaz, "stress Z", dim);
    }
};
}
#include "propertydefaults.hh"
#endif
