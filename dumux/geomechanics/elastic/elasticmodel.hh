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

#include "elasticproperties.hh"

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
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum {
        dim = GridView::dimension
    };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

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
        unsigned numScv = this->gridView_().size(dim);
        unsigned numElements = this->gridView_().size(0);

        ScalarField &ux = *writer.allocateManagedBuffer(numScv);
        ScalarField &uy = *writer.allocateManagedBuffer(numScv);
        ScalarField &uz = *writer.allocateManagedBuffer(numScv);
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

        FVElementGeometry fvGeometry;
        VolumeVariables volVars;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator eEndIt = this->gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt)
        {
            if(eIt->partitionType() == Dune::InteriorEntity)
            {
            int eIdx = this->problem_().model().elementMapper().map(*eIt);
            rank[eIdx] = this->gridView_().comm().rank();

            fvGeometry.update(this->gridView_(), *eIt);
            elemBcTypes.update(this->problem_(), *eIt, fvGeometry);

            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                int globalIdx = this->dofMapper().map(*eIt, scvIdx, dim);

                volVars.update(sol[globalIdx],
                               this->problem_(),
                               *eIt,
                               fvGeometry,
                               scvIdx,
                               false);

                ux[globalIdx] = volVars.displacement(0);
                if (dim >= 2)
                    uy[globalIdx] = volVars.displacement(1);
                if (dim >= 3)
                    uz[globalIdx] = volVars.displacement(2);
              };

            // In the box method, the stress is evaluated on the FE-Grid. However, to get an
            // average apparent stress for the cell, all contributing stresses have to be interpolated.
            DimMatrix stress;

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                            *eIt,
                            fvGeometry,
                            false /* isOldSol? */);

            // loop over the faces
            for (int fIdx = 0; fIdx < fvGeometry.numScvf; fIdx++)
            {
                stress = 0.0;
                //prepare the flux calculations (set up and prepare geometry, FE gradients)
                FluxVariables fluxVars(this->problem_(),
                              *eIt,
                              fvGeometry,
                              fIdx,
                              elemVolVars);

                stress = fluxVars.sigma();
                stress /= fvGeometry.numScvf;

                // Add up stresses for each cell.
                // Beware the sign convention applied here: compressive stresses are negative
                sigmax[eIdx] += stress[0];
                if (dim >= 2)
                {
                    sigmay[eIdx] += stress[1];
                }
                if (dim == 3)
                {
                    sigmaz[eIdx] += stress[2];
                }
            }
            }
        }


        writer.attachVertexData(ux, "ux");
        if (dim >= 2)
            writer.attachVertexData(uy, "uy");
        if (dim == 3)
            writer.attachVertexData(uz, "uz");
        writer.attachCellData(sigmax, "stress X", dim);
        if (dim >= 2)
        writer.attachCellData(sigmay, "stress Y", dim);
        if (dim == 3)
        writer.attachCellData(sigmaz, "stress Z", dim);
    }
};
}
#include "elasticpropertydefaults.hh"
#endif
