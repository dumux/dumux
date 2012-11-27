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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/

/*!
* \file
*
* \brief Adaption of the box scheme to the two-phase flow in discrete fracture-matrix 
*        model.
*/

#ifndef DUMUX_BOXMODELS_2PDFM_MODEL_HH
#define DUMUX_BOXMODELS_2PDFM_MODEL_HH

#include <dumux/boxmodels/2p/2pmodel.hh>

#include "2pdfmproperties.hh"
#include "2pdfmlocalresidual.hh"
#include "2pdfmproblem.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPDFMBoxModel
 * \brief A two-phase, isothermal flow model using the box scheme.
 *
 * This model implements two-phase flow of two immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multiphase Darcy
 * approach as the equation for the conservation of momentum, i.e.
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 \right\} - q_\alpha = 0 \;,
 \f]
 *
 * These equations are discretized by a fully-coupled vertex centered finite volume
 * (box) scheme as spatial and the implicit Euler method as time
 * discretization.
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
class TwoPDFMModel : public TwoPModel<TypeTag>
{
    typedef TwoPModel<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        pressureIdx = Indices::pressureIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::ctype CoordScalar;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param globalIdx The global index of the vertex
     * \param pvIdx The index of the primary variable
     */
    DUNE_DEPRECATED
    Scalar primaryVarWeight(int globalIdx, int pvIdx) const
    {
        if (pressureIdx == pvIdx)
        {
            return std::min(10.0 / this->prevSol()[globalIdx][pvIdx], 1.0);
        }
        
        return 1.0;
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     * \param sol The global solution vector
     * \param writer The writer for multi-file VTK datasets
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol, MultiWriter &writer)
    {
        bool velocityOutput = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddVelocity);
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

        // get the number of degrees of freedom and the number of vertices
        unsigned numDofs = this->numDofs();
        unsigned numVertices = this->problem_().gridView().size(dim);
    
        // velocity output currently only works for vertex data
        if (numDofs != numVertices)
        {
            velocityOutput = false;
        }
    
        // create the required scalar fields
        ScalarField *pW = writer.allocateManagedBuffer(numDofs);
        ScalarField *pN = writer.allocateManagedBuffer(numDofs);
        ScalarField *pC = writer.allocateManagedBuffer(numDofs);
        ScalarField *Sw = writer.allocateManagedBuffer(numDofs);
        ScalarField *Sn = writer.allocateManagedBuffer(numDofs);
        ScalarField *rhoW = writer.allocateManagedBuffer(numDofs);
        ScalarField *rhoN = writer.allocateManagedBuffer(numDofs);
        ScalarField *mobW = writer.allocateManagedBuffer(numDofs);
        ScalarField *mobN = writer.allocateManagedBuffer(numDofs);
        ScalarField *poro = writer.allocateManagedBuffer(numDofs);
        ScalarField *Te = writer.allocateManagedBuffer(numDofs);
        ScalarField *cellNum =writer.allocateManagedBuffer (numDofs);
        VectorField *velocityN = writer.template allocateManagedBuffer<double, dim>(numDofs);
        VectorField *velocityW = writer.template allocateManagedBuffer<double, dim>(numDofs);

        if(velocityOutput) // check if velocity output is demanded
        {
            // initialize velocity fields
            for (unsigned int i = 0; i < numVertices; ++i)
            {
                (*velocityN)[i] = Scalar(0);
                (*velocityW)[i] = Scalar(0);
                (*cellNum)[i] = Scalar(0.0);
            }
        }
        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer(numElements);

        FVElementGeometry fvGeometry;
        VolumeVariables volVars;
        ElementVolumeVariables elemVolVars;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            if(velocityOutput && !elemIt->geometry().type().isCube()){
                DUNE_THROW(Dune::InvalidStateException,
                           "Currently, velocity output only works for cubes. "
                           "Please set the EnableVelocityOutput property to false!");
            }
            int idx = this->elementMapper().map(*elemIt);
            (*rank)[idx] = this->gridView_().comm().rank();

            fvGeometry.update(this->gridView_(), *elemIt);

            if (numDofs == numElements) // element data
            {
                volVars.update(sol[idx],
                               this->problem_(),
                               *elemIt,
                               fvGeometry,
                               /*scvIdx=*/0,
                               false);

                (*pW)[idx] = volVars.pressure(wPhaseIdx);
                (*pN)[idx] = volVars.pressure(nPhaseIdx);
                (*pC)[idx] = volVars.capillaryPressure();
                (*Sw)[idx] = volVars.saturation(wPhaseIdx);
                (*Sn)[idx] = volVars.saturation(nPhaseIdx);
                (*rhoW)[idx] = volVars.density(wPhaseIdx);
                (*rhoN)[idx] = volVars.density(nPhaseIdx);
                (*mobW)[idx] = volVars.mobility(wPhaseIdx);
                (*mobN)[idx] = volVars.mobility(nPhaseIdx);
                (*poro)[idx] = volVars.porosity();
                (*Te)[idx] = volVars.temperature();
            }
            else // vertex data
            {
                int numVerts = elemIt->template count<dim> ();
                for (int scvIdx = 0; scvIdx < numVerts; ++scvIdx)
                {
                    int globalIdx = this->vertexMapper().map(*elemIt, scvIdx, dim);
                    volVars.update(sol[globalIdx],
                                   this->problem_(),
                                   *elemIt,
                                   fvGeometry,
                                   scvIdx,
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
                    if(velocityOutput)
                    {
                        (*cellNum)[globalIdx] += 1;
                    }
                }

                if(velocityOutput)
                {
                    // calculate vertex velocities
                    GlobalPosition tmpVelocity[numPhases];

                    for(int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    {
                        tmpVelocity[phaseIdx]  = Scalar(0.0);
                    }

                    typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > SCVVelocities;
                    SCVVelocities scvVelocityW(8), scvVelocityN(8);

                    scvVelocityW = 0;
                    scvVelocityN = 0;

                    ElementVolumeVariables elemVolVars;

                    elemVolVars.update(this->problem_(),
                                       *elemIt,
                                       fvGeometry,
                                       false /* oldSol? */);

                    for (int faceIdx = 0; faceIdx < fvGeometry.numEdges; faceIdx++)
                    {
                        FluxVariables fluxVars(this->problem_(),
                                             *elemIt,
                                             fvGeometry,
                                             faceIdx,
                                             elemVolVars);

                        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                        {
                            // local position of integration point
                            const Dune::FieldVector<Scalar, dim>& localPosIP = fvGeometry.subContVolFace[faceIdx].ipLocal;

                            // Transformation of the global normal vector to normal vector in the reference element
                            const typename Element::Geometry::JacobianTransposed& jacobianT1 = 
                                elemIt->geometry().jacobianTransposed(localPosIP);

                            const GlobalPosition globalNormal = fluxVars.face().normal;
                            GlobalPosition localNormal;
                            jacobianT1.mv(globalNormal, localNormal);
                            // note only works for cubes
                            const Scalar localArea = pow(2,-(dim-1));

                            localNormal /= localNormal.two_norm();

                            // Get the Darcy velocities. The Darcy velocities are divided by the area of the subcontrolvolumeface
                            // in the reference element.
                            PhasesVector q;
                            q[phaseIdx] = fluxVars.volumeFlux(phaseIdx) / localArea;

                            // transform the normal Darcy velocity into a vector
                            tmpVelocity[phaseIdx] = localNormal;
                            tmpVelocity[phaseIdx] *= q[phaseIdx];

                            if(phaseIdx == wPhaseIdx)
                            {
                                scvVelocityW[fluxVars.face().i] += tmpVelocity[phaseIdx];
                                scvVelocityW[fluxVars.face().j] += tmpVelocity[phaseIdx];
                            }
                            else if(phaseIdx == nPhaseIdx)
                            {
                                scvVelocityN[fluxVars.face().i] += tmpVelocity[phaseIdx];
                                scvVelocityN[fluxVars.face().j] += tmpVelocity[phaseIdx];
                            }
                        }
                    }
                    typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
                    const Dune::FieldVector<Scalar, dim> &localPos
                        = ReferenceElements::general(elemIt->geometry().type()).position(0, 0);

                    // get the transposed Jacobian of the element mapping
                    const typename Element::Geometry::JacobianTransposed& jacobianT2
                        = elemIt->geometry().jacobianTransposed(localPos);

                    // transform vertex velocities from local to global coordinates
                    for (int i = 0; i < numVerts; ++i)
                    {
                        int globalIdx = this->vertexMapper().map(*elemIt, i, dim);
                        // calculate the subcontrolvolume velocity by the Piola transformation
                        Dune::FieldVector<CoordScalar, dim> scvVelocity(0);

                        jacobianT2.mtv(scvVelocityW[i], scvVelocity);
                        scvVelocity /= elemIt->geometry().integrationElement(localPos);
                        // add up the wetting phase subcontrolvolume velocities for each vertex
                        (*velocityW)[globalIdx] += scvVelocity;

                        jacobianT2.mtv(scvVelocityN[i], scvVelocity);
                        scvVelocity /= elemIt->geometry().integrationElement(localPos);
                        // add up the nonwetting phase subcontrolvolume velocities for each vertex
                        (*velocityN)[globalIdx] += scvVelocity;
                    }
                }
            }
        }

        if (numDofs == numElements) // element data
        {
            writer.attachCellData(*Sn, "Sn");
            writer.attachCellData(*Sw, "Sw");
            writer.attachCellData(*pN, "pn");
            writer.attachCellData(*pW, "pw");
            writer.attachCellData(*pC, "pc");
            writer.attachCellData(*rhoW, "rhoW");
            writer.attachCellData(*rhoN, "rhoN");
            writer.attachCellData(*mobW, "mobW");
            writer.attachCellData(*mobN, "mobN");
            writer.attachCellData(*poro, "porosity");
            writer.attachCellData(*Te, "temperature");
        }
        else // vertex data
        {
            writer.attachVertexData(*Sn, "Sn");
            writer.attachVertexData(*Sw, "Sw");
            writer.attachVertexData(*pN, "pn");
            writer.attachVertexData(*pW, "pw");
            writer.attachVertexData(*pC, "pc");
            writer.attachVertexData(*rhoW, "rhoW");
            writer.attachVertexData(*rhoN, "rhoN");
            writer.attachVertexData(*mobW, "mobW");
            writer.attachVertexData(*mobN, "mobN");
            writer.attachVertexData(*poro, "porosity");
            writer.attachVertexData(*Te, "temperature");
            if(velocityOutput) // check if velocity output is demanded
            {
                // divide the vertex velocities by the number of adjacent scvs i.e. cells
                for(unsigned int globalIdx = 0; globalIdx < numVertices; ++globalIdx)
                {
                    (*velocityW)[globalIdx] /= (*cellNum)[globalIdx];
                    (*velocityN)[globalIdx] /= (*cellNum)[globalIdx];
                }
                writer.attachVertexData(*velocityW,  "velocityW", dim);
                writer.attachVertexData(*velocityN,  "velocityN", dim);
            }
        }
        writer.attachCellData(*rank, "process rank");
    }
};
} // end namespace

#include "2pdfmpropertydefaults.hh"

#endif // DUMUX_BOXMODELS_2PDFM_MODEL_HH
