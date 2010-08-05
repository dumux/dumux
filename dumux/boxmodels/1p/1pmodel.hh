// $Id: 1pboxmodel.hh 3738 2010-06-15 14:01:09Z lauser $
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_1P_MODEL_HH
#define DUMUX_1P_MODEL_HH

#include <dumux/boxmodels/common/boxmodel.hh>

#include "1plocalresidual.hh"
#include "1pproblem.hh"

namespace Dumux
{
/*!
 * \ingroup BoxProblems
 * \defgroup OnePBoxProblems One-phase box problems
 */

/*!
 * \ingroup BoxModels
 * \defgroup OnePBoxModel One-phase box model
 */

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
 * single phase flow modeling, if in the problem file a fluid with constant density is chosen.
 */
template<class TypeTag >
class OnePBoxModel : public BoxModel<TypeTag>
{
    typedef OnePBoxModel<TypeTag> ThisType;
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementMapper)) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePIndices)) Indices;

    enum {
        dim = GridView::dimension,
    };

public:
    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_().gridView().size(dim);
        ScalarField *p = writer.template createField<Scalar, 1> (numVertices);

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank =
                writer.template createField<Scalar, 1> (numElements);

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

                (*p)[globalIdx] = volVars.pressure;
            };
        }

        writer.addVertexData(p, "p");
        writer.addCellData(rank, "process rank");
    }
};
}

#endif
