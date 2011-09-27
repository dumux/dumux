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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;


    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),

        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx,
    };
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::ctype CoordScalar;

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
        bool velocityOutput = GET_PROP_VALUE(TypeTag, PTAG(EnableVelocityOutput));
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

        // create the required scalar fields
        unsigned numVertices = this->problem_().gridView().size(dim);
        ScalarField *pW = writer.allocateManagedBuffer(numVertices);
        ScalarField *pN = writer.allocateManagedBuffer(numVertices);
        ScalarField *pC = writer.allocateManagedBuffer(numVertices);
        ScalarField *Sw = writer.allocateManagedBuffer(numVertices);
        ScalarField *Sn = writer.allocateManagedBuffer(numVertices);
        ScalarField *rhoW = writer.allocateManagedBuffer(numVertices);
        ScalarField *rhoN = writer.allocateManagedBuffer(numVertices);
        ScalarField *mobW = writer.allocateManagedBuffer(numVertices);
        ScalarField *mobN = writer.allocateManagedBuffer(numVertices);
        ScalarField *poro = writer.allocateManagedBuffer(numVertices);
        ScalarField *Te = writer.allocateManagedBuffer(numVertices);
        ScalarField *cellNum =writer.allocateManagedBuffer (numVertices);
        VectorField *velocityN = writer.template allocateManagedBuffer<double, dim>(numVertices);
        VectorField *velocityW = writer.template allocateManagedBuffer<double, dim>(numVertices);

        if(velocityOutput) // check if velocity output is demanded
        {
            // initialize velocity fields
            for (int i = 0; i < numVertices; ++i)
            {
            	(*velocityN)[i] = Scalar(0);
            	(*velocityW)[i] = Scalar(0);
            	(*cellNum)[i] = Scalar(0.0);
            }
        }
        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer(numElements);

        FVElementGeometry fvElemGeom;
        VolumeVariables volVars;
        ElementVolumeVariables elemVolVars;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
#warning "currently, velocity output only works for cubes and is set to false for simplices"
        	if(elemIt->geometry().type().isCube() == false){
        		velocityOutput = false;
        	}
            int idx = this->elementMapper().map(*elemIt);
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
                if(velocityOutput)
                {
                    (*cellNum)[globalIdx] += 1;
                }
            };
            
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
                    	          fvElemGeom,
                        	      false /* oldSol? */);
        		
                for (int faceIdx = 0; faceIdx< fvElemGeom.numEdges; faceIdx++)
                {

                    FluxVariables fluxDat(this->problem_(),
                                  *elemIt,
                                  fvElemGeom,
                                  faceIdx,
                                  elemVolVars);

                    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    {

     					// data attached to upstream and the downstream vertices
                        // of the current phase
                        const VolumeVariables up =
                            elemVolVars[fluxDat.upstreamIdx(phaseIdx)];
                        const VolumeVariables dn =
                            elemVolVars[fluxDat.downstreamIdx(phaseIdx)];

                      // calculate the flux in the normal direction of the
                      // current sub control volume face
                      GlobalPosition tmpVec(0);
                      fluxDat.intrinsicPermeability().mv(fluxDat.potentialGrad(phaseIdx),
                                                     tmpVec);
                      const GlobalPosition globalNormal = fluxDat.face().normal;
                      const Scalar normalFlux = - (tmpVec*globalNormal);

                      // local position of integration point
                      const Dune::FieldVector<Scalar, dim>& localPosIP = fvElemGeom.subContVolFace[faceIdx].ipLocal;

                      // Transformation of the global normal vector to normal vector in the reference element
                      const Dune::FieldMatrix<CoordScalar, dim, dim> jacobianT1 = elemIt->geometry().jacobianTransposed(localPosIP);

                      GlobalPosition localNormal(0);
                      jacobianT1.mv(globalNormal, localNormal);
                  	  // note only works for cubes
                      const Scalar localArea = pow(2,-(dim-1));

                      localNormal /= localNormal.two_norm();

                      // Get the Darcy velocities. The Darcy velocities are divided by the area of the subcontrolvolumeface
                      // in the reference element.
                      const Scalar massUpwindWeight = GET_PARAM(TypeTag, Scalar, MassUpwindWeight);
                      PhasesVector q;
                      q[phaseIdx] = normalFlux
                                		   * (massUpwindWeight
                                		   * up.mobility(phaseIdx)
                                		   + (1- massUpwindWeight)
                                		   * dn.mobility(phaseIdx)) / localArea;

                      // transform the normal Darcy velocity into a vector
                      tmpVelocity[phaseIdx] = localNormal;
                      tmpVelocity[phaseIdx] *= q[phaseIdx];

                      if(phaseIdx == wPhaseIdx){
                      scvVelocityW[fluxDat.face().i] += tmpVelocity[phaseIdx];
                      scvVelocityW[fluxDat.face().j] += tmpVelocity[phaseIdx];
                      }
                      else if(phaseIdx == nPhaseIdx){
                      scvVelocityN[fluxDat.face().i] += tmpVelocity[phaseIdx];
                      scvVelocityN[fluxDat.face().j] += tmpVelocity[phaseIdx];
                      }
                   }
                }
                typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
                const Dune::FieldVector<Scalar, dim> &localPos 
                    = ReferenceElements::general(elemIt->geometry().type()).position(0, 0);

     			// get the transposed Jacobian of the element mapping
                const Dune::FieldMatrix<CoordScalar, dim, dim> &jacobianT2
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
            if(velocityOutput)
            {
            	// divide the vertex velocities by the number of adjacent scvs i.e. cells
    	        for(int globalIdx = 0; globalIdx<numVertices; ++globalIdx){
    	        (*velocityW)[globalIdx] /= (*cellNum)[globalIdx];
    	        (*velocityN)[globalIdx] /= (*cellNum)[globalIdx];
    	        }
            }

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
        	writer.attachVertexData(*velocityW,  "velocityW", dim);
        	writer.attachVertexData(*velocityN,  "velocityN", dim);
        }
        writer.attachCellData(*rank, "process rank");
    }
};
}

#include "2ppropertydefaults.hh"

#endif
