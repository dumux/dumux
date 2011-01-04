// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
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
 * \brief Base class for all models which use the single-phase,
 *        two-component box model
 *        Adaption of the BOX scheme to the one-phase two-component flow model.
 */


#ifndef DUMUX_ONEP_TWOC_MODEL_HH
#define DUMUX_ONEP_TWOC_MODEL_HH

#include "1p2cproperties.hh"
#include "1p2cproblem.hh"

#include <dumux/boxmodels/common/boxmodel.hh>

namespace Dumux
{

/*!
 * \ingroup BoxModels
 * \defgroup OnePTwoCBoxModel One-phase Two-component box model
 */

/*!
 * \ingroup OnePTwoCBoxModel
 * \brief Adaption of the BOX scheme to the one-phase two-component flow model.
 *
 * This model implements an one-phase flow of an incompressible fluid, that consists of two components,
 * using a standard Darcy
 * approach as the equation for the conservation of momentum:
 \f[
 v_{D} = - \frac{K}{\mu}
 \left(\text{grad} p - \varrho g \right)
 \f]
 *
 * Gravity can be enabled or disabled via the Property system.
 * By inserting this into the continuity equation, one gets
 \f[
 - \text{div} \left\{
   \varrho \frac{K}{\mu}  \left(\text{grad} p - \varrho g \right)
 \right\} = q \;,
 \f]
 *
 * The transport of the components is described by the following equation:
 \f[
 \Phi \varrho \frac{ \partial x}{\partial t} - \text{div} \left( \varrho \frac{K x}{\mu} \left( \text{grad} p -
 \varrho g \right) + \varrho \tau \Phi D \text{grad} x \right) = q.
 \f]
 *
 * All equations are discretized using a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as time discretization.
 *
 * The primary variables are the pressure \f$p\f$ and the mole fraction of dissolved component \f$x\f$.
 */

template<class TypeTag >
class OnePTwoCBoxModel : public BoxModel<TypeTag>
{
    typedef OnePTwoCBoxModel<TypeTag> ThisType;
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
     * \brief \copybrief Dumux::BoxModel::addOutputVtkFields
     *
     * \copydetails Dumux::BoxModel::addOutputVtkFields
     *
     * Specialization for the OnePTwoCBoxModel, adding pressure,
     * mass and mole fractions, and the process rank to the VTK writer.
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_().gridView().size(dim);
        ScalarField *pressure = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *moleFrac0 = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *moleFrac1 = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massFrac0 = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massFrac1 = writer.template createField<Scalar, 1>(numVertices);

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


                (*pressure)[globalIdx] = volVars.pressure();
                (*moleFrac0)[globalIdx] = volVars.moleFrac(0);
                (*moleFrac1)[globalIdx] = volVars.moleFrac(1);
                (*massFrac0)[globalIdx] = volVars.massFrac(0);
                (*massFrac1)[globalIdx] = volVars.massFrac(1);
            };
        }

        writer.addVertexData(pressure, "p");
        writer.addVertexData(moleFrac0, "x_0");
        writer.addVertexData(moleFrac1, "x_1");
        writer.addVertexData(massFrac0, "X_0");
        writer.addVertexData(massFrac1, "X_1");
        writer.addCellData(rank, "process rank");
    }

};
}

#include "1p2cpropertydefaults.hh"

#endif
