// $Id$
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2007-2007 by Bernd Flemisch                               *
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
* \brief Adaption of the box scheme to the two-phase flow model.
*/

#ifndef DUMUX_TWOP_MODEL_HH
#define DUMUX_TWOP_MODEL_HH

#include "2pproperties.hh"
#include "2plocalresidual.hh"
#include "2pproblem.hh"

namespace Dumux
{

/*!
 * \ingroup BoxModels
 * \defgroup TwoPBoxModel Two-phase box model
 */

/*!
 * \ingroup TwoPBoxModel
 * \brief Adaption of the BOX scheme to the twophase flow model.
 *
 * This model implements two-phase flow of two completely immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multiphase Darcy
 * approach as the equation for the conservation of momentum:
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad} p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K} \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 \right\} - q_\alpha = 0 \;,
 \f]
 * discretized by a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$, the number of
 * unknowns can be reduced to two. Currently the model supports
 * choosing either \f$p_w\f$ and \f$S_n\f$ or \f$p_n\f$ and \f$S_w\f$
 * as primary variables. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * <tt>TwoPCommonIndices::pWsN</tt> or <tt>TwoPCommonIndices::pNsW</tt>. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 */
template<class TypeTag >
class TwoPModel : public BoxModel<TypeTag>
{
    typedef TwoPModel<TypeTag> ThisType;
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    enum {
        dim = GridView::dimension,
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx,
    };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

public:
    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     *        \param globalVertexIdx The global index of the vertex
     *        \param pvIdx The index of the primary variable
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        if (Indices::pressureIdx == pvIdx)
            return std::min(10.0/this->prevSol()[globalVertexIdx][pvIdx], 1.0);
        return 1;
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     *        \param sol The global solution vector
     *        \param writer The writer for multi-file VTK datasets
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_().gridView().size(dim);
        ScalarField *pW = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *pN = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *pC = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *Sw = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *Sn = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *rhoW = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *rhoN = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *mobW = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *mobN = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *poro = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *Te = writer.template createField<Scalar, 1> (numVertices);

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank =
                writer.template createField<Scalar, 1> (numElements);

        FVElementGeometry fvElemGeom;
        VolumeVariables volVars;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            int idx = this->problem_().model().elementMapper().map(*elemIt);
            (*rank)[idx] = this->gridView_().comm().rank();

            fvElemGeom.update(this->gridView_(), *elemIt);

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

                (*pW)[globalIdx] = volVars.pressure(wPhaseIdx);
                (*pN)[globalIdx] = volVars.pressure(nPhaseIdx);
                (*pC)[globalIdx] = volVars.capillaryPressure();
                (*Sw)[globalIdx] = volVars.saturation(wPhaseIdx);
                (*Sn)[globalIdx] = volVars.saturation(nPhaseIdx);
                (*rhoW)[globalIdx] = volVars.density(wPhaseIdx);
                (*rhoN)[globalIdx] = volVars.density(nPhaseIdx);
                (*mobW)[globalIdx] = volVars.mobility(wPhaseIdx);
                (*mobN)[globalIdx] = volVars.mobility(nPhaseIdx);
                (*poro)[globalIdx] = volVars.porosity();
                (*Te)[globalIdx] = volVars.temperature();
            };
        }

        writer.addVertexData(Sn, "Sn");
        writer.addVertexData(Sw, "Sw");
        writer.addVertexData(pN, "pn");
        writer.addVertexData(pW, "pw");
        writer.addVertexData(pC, "pc");
        writer.addVertexData(rhoW, "rhoW");
        writer.addVertexData(rhoN, "rhoN");
        writer.addVertexData(mobW, "mobW");
        writer.addVertexData(mobN, "mobN");
        writer.addVertexData(poro, "porosity");
        writer.addVertexData(Te, "temperature");
        writer.addCellData(rank, "process rank");
    }
};
}

#include "2ppropertydefaults.hh"

#endif
