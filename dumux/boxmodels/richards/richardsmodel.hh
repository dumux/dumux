// $Id$
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
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
#ifndef DUMUX_RICHARDS_MODEL_HH
#define DUMUX_RICHARDS_MODEL_HH

#include <dumux/boxmodels/common/boxmodel.hh>

#include "richardslocalresidual.hh"
#include "richardsproblem.hh"

namespace Dumux
{
/*!
 * \ingroup BoxProblems
 * \defgroup RichardsBoxProblems Richards box problems
 */

/*!
 * \ingroup BoxModels
 * \defgroup RichardsModel Richards box model
 */

/*!
 * \ingroup RichardsModel
 * \brief Adaption of the BOX scheme to the isothermal Richards model.
 *
 *
 * In the unsaturated zone, Richards' equation can be used.
 * Gas has resistance against the water flow in porous media.
 * However, viscosity of air is about 1\% of the viscosity of water,
 * which makes it highly mobile compared to the water phase.
 * Therefore, in Richards` equation only water phase with capillary effects are considered,
 * where pressure of the gas phase is set to a reference pressure (\f${p_n}_{ref}\f$).
 *
 * \f{align*}
 * \varrho \hspace{1mm} \phi \hspace{1mm} \frac{\partial S_w}{\partial p_c} \frac{\partial p_c}{\partial t} - \nabla \cdot (\frac{kr_w}{\mu_w} \hspace{1mm} \varrho_w \hspace{1mm} K \hspace{1mm}
 * (\nabla p_w - \varrho_w \hspace{1mm} \vec{g})) \hspace{1mm} = \hspace{1mm} q,
 * \f}
 * where \f$p_w = {p_n}_{ref} - p_c\f$.
 * Here \f$ p_w \f$, \f$ p_c \f$, and \f$ {p_n}_{ref} \f$
 * denote water pressure, capillary pressure, and non-wetting phase reference pressure, repectively.
 *
 * To overcome convergence problems, \f$ \frac{\partial S_w}{\partial p_c} \f$ is taken from the old iteration step.
 *
 */
template<class TypeTag >
class RichardsModel : public BoxModel<TypeTag>
{
    typedef RichardsModel<TypeTag> ThisType;
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

    enum {
        dim = GridView::dimension,
    };

public:
    /*!
     * \brief All relevant primary and secondary of a given
     *        solution to an ouput writer.
     */
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol, MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);
        ScalarField *Sw = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sn = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pC = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pW = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *dSwdpC = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *rhoW = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *mobW = writer.template createField<Scalar, 1>(numVertices);

        FVElementGeometry fvElemGeom;
        VolumeVariables volVars;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
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

                (*Sw)[globalIdx] = volVars.Sw;
                (*Sn)[globalIdx] = 1.0 - volVars.Sw;
                (*pC)[globalIdx] = volVars.pC;
                (*pW)[globalIdx] = volVars.pW;
                (*dSwdpC)[globalIdx] = volVars.dSwdpC;
                (*rhoW)[globalIdx] = volVars.densityW;
                (*mobW)[globalIdx] = volVars.mobilityW;
            };
        }


        writer.addVertexData(Sw, "Sw");
        writer.addVertexData(Sn, "Sn");
        writer.addVertexData(pC, "pC");
        writer.addVertexData(pW, "pW");
        writer.addVertexData(dSwdpC, "dSwdpC");
        writer.addVertexData(rhoW, "rhoW");
        writer.addVertexData(mobW, "mobW");
    }
};
}

#endif
