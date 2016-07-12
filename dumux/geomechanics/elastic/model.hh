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
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dim = GridView::dimension };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
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
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > VectorField;

        // create the required scalar fields
        unsigned numDofs = this->numDofs();
        unsigned numElements = this->gridView_().size(0);

        ScalarField &ux = *writer.allocateManagedBuffer(numDofs);
        ScalarField &uy = *writer.allocateManagedBuffer(numDofs);
        ScalarField &uz = *writer.allocateManagedBuffer(numDofs);
        VectorField &sigmax = *writer.template allocateManagedBuffer<Scalar, dim>(numElements);
        VectorField &sigmay = *writer.template allocateManagedBuffer<Scalar, dim>(numElements);
        VectorField &sigmaz = *writer.template allocateManagedBuffer<Scalar, dim>(numElements);

        // initialize stress fields
        for (unsigned int i = 0; i < numElements; ++i)
        {
            sigmax[i] = 0;
            if (dim > 1)
            {
                sigmay[i] = 0;
            }
            if (dim > 2)
            {
                sigmaz[i] = 0;
            }
         }

        ScalarField &rank = *writer.allocateManagedBuffer(numElements);

        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
        {
            int eIdx = this->problem_().model().elementMapper().index(element);
            rank[eIdx] = this->gridView_().comm().rank();

            // make sure FVElementGeometry and the volume variables are bound to the element
            this->fvGeometries_().bind(element);
            this->curVolVars_().bind(element);

            const auto& fvGeometry = this->fvGeometries(element);
            for (const auto& scv : fvGeometry.scvs())
            {
                int dofIdxGlobal = scv.dofIndex();
                const auto& volVars = this->curVolVars(scv);

                ux[dofIdxGlobal] = volVars.displacement(0);
                if (dim >= 2)
                    uy[dofIdxGlobal] = volVars.displacement(1);
                if (dim >= 3)
                    uz[dofIdxGlobal] = volVars.displacement(2);
            }

            // In the box method, the stress is evaluated on the FE-Grid. However, to get an
            // average apparent stress for the cell, all contributing stresses have to be interpolated.
            DimMatrix stress(0.0);
            unsigned int counter = 0;

            // loop over the faces
            for (const auto& scvFace : fvGeometry.scvfs())
            {
                if (scvFace.boundary())
                {
                    BoundaryTypes bcTypes = this->problem_().boundaryTypes(element, scvFace);
                    if (bcTypes.hasNeumann())
                        continue;
                }

                //prepare the flux calculations (set up and prepare geometry, FE gradients)
                FluxVariables fluxVars;
                fluxVars.initAndComputeFluxes(this->problem_(), element, scvFace);

                // Add up stresses for each scv face.
                // Beware the sign convention applied here: compressive stresses are negative
                stress += fluxVars.stressTensor();
                counter++;
            }

            // divide by the number of added stress tensors and add to container
            stress /= counter;
            sigmax[eIdx] += stress[0];
            if (dim >= 2)
                sigmay[eIdx] += stress[1];
            if (dim == 3)
                sigmaz[eIdx] += stress[2];
            }
        }

        writer.attachDofData(ux, "ux", isBox);
        if (dim >= 2)
            writer.attachDofData(uy, "uy", isBox);
        if (dim == 3)
            writer.attachDofData(uz, "uz", isBox);
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
