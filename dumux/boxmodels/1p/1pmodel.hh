// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Base class for all models which use the one-phase,
 *        box model.
 *        Adaption of the BOX scheme to the one-phase flow model.
 */

#ifndef DUMUX_1P_MODEL_HH
#define DUMUX_1P_MODEL_HH

#include <dumux/boxmodels/common/boxmodel.hh>

#include "1plocalresidual.hh"
#include "1pproblem.hh"

namespace Dumux
{
/*!
 * \ingroup OnePBoxModel
 * \brief Adaption of the BOX scheme to the single phase isothermal flow model.
 *
 * Single phase compressible isothermal flow model,
 * \f{align*}
 * \phi \frac{\partial \varrho}{\partial t} + \vec{\nabla} \cdot (- \varrho \frac{\bar{\bar{K}}}{\mu} ( \nabla p -\varrho \vec{g})) = q,
 * \f}
 * discretized using a vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as time discretization.
 * Of course, the model can also be used for incompressible
 * single phase flow modeling, if a fluid with constant density is chosen in the problem file.
 */
template<class TypeTag >
class OnePBoxModel : public BoxModel<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    enum { dim = GridView::dimension };

public:
    /*!
     * \brief \copybrief Dumux::BoxModel::addOutputVtkFields
     *
     * \copydetails Dumux::BoxModel::addOutputVtkFields
     *
     * Specialization for the OnePBoxModel, adding the pressure and
     * the process rank to the VTK writer.
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_().gridView().size(dim);
        ScalarField *p = writer.allocateManagedBuffer(numVertices);
        ScalarField *K = writer.allocateManagedBuffer(numVertices);

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer(numElements);

        FVElementGeometry fvElemGeom;
        VolumeVariables volVars;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            int idx = this->problem_().model().elementMapper().map(*elemIt);
            (*rank)[idx] = this->gridView_().comm().rank();

            fvElemGeom.update(this->gridView_(), *elemIt);
            elemBcTypes.update(this->problem_(), *elemIt, fvElemGeom);

            int numVerts = elemIt->template count<dim> ();
            for (int i = 0; i < numVerts; ++i)
            {
                int globalIdx = this->vertexMapper().map(*elemIt, i, dim);
                volVars.update(sol[globalIdx],
                               this->problem_(),
                               *elemIt,
                               fvElemGeom,
                               i,
                               false);
                const SpatialParameters &spatialParams = this->problem_().spatialParameters();

                (*p)[globalIdx] = volVars.pressure();
                (*K)[globalIdx] = spatialParams.intrinsicPermeability(*elemIt,
                                                                    fvElemGeom,
                                                                    i);
            };
        }

        writer.attachVertexData(*p, "p");
        writer.attachVertexData(*K, "K");
        writer.attachCellData(*rank, "process rank");
    }
};
}

#include "1ppropertydefaults.hh"

#endif
