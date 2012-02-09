// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief Base class for all models which use the single-phase,
 *        two-component box model.
 *        Adaption of the BOX scheme to the one-phase two-component flow model.
 */

#ifndef DUMUX_ONEP_TWOC_MODEL_HH
#define DUMUX_ONEP_TWOC_MODEL_HH

#include "1p2cproperties.hh"
#include "1p2cproblem.hh"
#include "1p2clocalresidual.hh"

#include <dumux/boxmodels/common/boxmodel.hh>

namespace Dumux
{

/*!
 * \ingroup OnePTwoCBoxModel
 * \brief Adaption of the BOX scheme to the one-phase two-component flow model.
 *
 * This model implements a one-phase flow of a compressible fluid, that consists of two components,
 * using a standard Darcy
 * approach as the equation for the conservation of momentum:
 \f[
 v_{D} = - \frac{\textbf K}{\mu}
 \left(\text{grad} p - \varrho {\textbf g} \right)
 \f]
 *
 * Gravity can be enabled or disabled via the property system.
 * By inserting this into the continuity equation, one gets
 \f[
 \Phi \frac{\partial \varrho}{\partial t} - \text{div} \left\{
   \varrho \frac{\textbf K}{\mu}  \left(\text{grad}\, p - \varrho {\textbf g} \right)
 \right\} = q \;,
 \f]
 *
 * The transport of the components is described by the following equation:
 \f[
 \Phi \frac{ \partial \varrho x}{\partial t} - \text{div} \left( \varrho \frac{{\textbf K} x}{\mu} \left( \text{grad}\, p -
 \varrho {\textbf g} \right) + \varrho \tau \Phi D \text{grad} x \right) = q.
 \f]
 *
 * All equations are discretized using a fully-coupled vertex-centered
 * finite volume (box) scheme as spatial and
 * the implicit Euler method as time discretization.
 *
 * The primary variables are the pressure \f$p\f$ and the mole or mass fraction of dissolved component \f$x\f$.
 */

template<class TypeTag >
class OnePTwoCBoxModel : public BoxModel<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    OnePTwoCBoxModel()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        upwindWeight_ = GET_PARAM(TypeTag, Scalar, UpwindWeight);
    }

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
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_().gridView().size(dim);
        ScalarField &pressure = *writer.allocateManagedBuffer(numVertices);
        ScalarField &delp = *writer.allocateManagedBuffer(numVertices);
        ScalarField &moleFrac0 = *writer.allocateManagedBuffer(numVertices);
        ScalarField &moleFrac1 = *writer.allocateManagedBuffer(numVertices);
        ScalarField &massFrac0 = *writer.allocateManagedBuffer(numVertices);
        ScalarField &massFrac1 = *writer.allocateManagedBuffer(numVertices);
        ScalarField &rho = *writer.allocateManagedBuffer(numVertices);
        ScalarField &mu = *writer.allocateManagedBuffer(numVertices);
        ScalarField &delFrac= *writer.allocateManagedBuffer(numVertices);
#ifdef VELOCITY_OUTPUT // check if velocity output is demanded
        ScalarField &velocityX = *writer.allocateManagedBuffer(numVertices);
        ScalarField &velocityY = *writer.allocateManagedBuffer(numVertices);
        ScalarField &velocityZ = *writer.allocateManagedBuffer(numVertices);
        //use vertiacl faces for vx and horizontal faces for vy calculation
        std::vector<GlobalPosition> boxSurface(numVertices);
        // initialize velocity fields
          for (int i = 0; i < numVertices; ++i)
          {

              velocityX[i] = 0;
              if (dim > 1)
              {
                  velocityY[i] = 0;
              }
              if (dim > 2)
              {
                  velocityZ[i] = 0;
              }
              boxSurface[i] = Scalar(0.0); // initialize the boundary surface of the fv-boxes
          }
#endif
        unsigned numElements = this->gridView_().size(0);
        ScalarField &rank =
                *writer.allocateManagedBuffer(numElements);

        FVElementGeometry fvElemGeom;
        VolumeVariables volVars;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            int idx = this->problem_().model().elementMapper().map(*elemIt);
            rank[idx] = this->gridView_().comm().rank();

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

                pressure[globalIdx] = volVars.pressure();
                delp[globalIdx] = volVars.pressure() - 1e5;
                moleFrac0[globalIdx] = volVars.moleFraction(0);
                moleFrac1[globalIdx] = volVars.moleFraction(1);
                massFrac0[globalIdx] = volVars.massFraction(0);
                massFrac1[globalIdx] = volVars.massFraction(1);
                rho[globalIdx] = volVars.density();
                mu[globalIdx] = volVars.viscosity();
                delFrac[globalIdx] = volVars.massFraction(1)-volVars.moleFraction(1);
            };

#ifdef VELOCITY_OUTPUT // check if velocity output is demanded
            // In the box method, the velocity is evaluated on the FE-Grid. However, to get an
            // average apparent velocity at the vertex, all contributing velocities have to be interpolated.
            GlobalPosition velocity;

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               *elemIt,
                               fvElemGeom,
                               false /* isOldSol? */);
            // loop over the phases
            for (int faceIdx = 0; faceIdx < fvElemGeom.numEdges; faceIdx++)
            {
                velocity = 0.0;
                //prepare the flux calculations (set up and prepare geometry, FE gradients)
                FluxVariables fluxVars(this->problem_(),
                                       *elemIt,
                                       fvElemGeom,
                                       faceIdx,
                                       elemVolVars);

                //use vertiacl faces for vx and horizontal faces for vy calculation
                GlobalPosition xVector(0), yVector(0);
                xVector[0] = 1; yVector[1] = 1;

                Dune::SeqScalarProduct<GlobalPosition> sp;

                Scalar xDir = std::abs(sp.dot(fluxVars.face().normal, xVector));
                Scalar yDir = std::abs(sp.dot(fluxVars.face().normal, yVector));

                // up+downstream node
                const VolumeVariables &up =
                    elemVolVars[fluxVars.upstreamIdx()];
                const VolumeVariables &dn =
                    elemVolVars[fluxVars.downstreamIdx()];

                //get surface area to weight velocity at the IP with the surface area
                Scalar scvfArea = fluxVars.face().normal.two_norm();

                int vertIIdx = this->problem_().vertexMapper().map(
                    *elemIt, fluxVars.face().i, dim);
                int vertJIdx = this->problem_().vertexMapper().map(
                    *elemIt, fluxVars.face().j, dim);

                //use vertical faces (horizontal noraml vector) to calculate vx
                //in case of heterogeneities it seams to be better to define intrinisc permeability elementwise
                if(xDir > yDir)//(fluxVars.face().normal[0] > 1e-10 || fluxVars.face().normal[0] < -1e-10)// (xDir > yDir)
                {
                    // get darcy velocity
                    //calculate (v n) n/A
                    Scalar tmp = fluxVars.KmvpNormal();
                    velocity = fluxVars.face().normal;
                    velocity *= tmp;
                    velocity /= scvfArea;
                    velocity *= (upwindWeight_ / up.viscosity() +
                                 (1 - upwindWeight_)/ dn.viscosity());

                    // add surface area for weighting purposes
                    boxSurface[vertIIdx][0] += scvfArea;
                    boxSurface[vertJIdx][0] += scvfArea;

                    velocityX[vertJIdx] += velocity[0];
                    velocityX[vertIIdx] += velocity[0];

                }
                if (yDir > xDir)//(fluxVars.face().normal[1] > 1e-10 || fluxVars.face().normal[1] < -1e-10)// (yDir > xDir)
                {
                    // get darcy velocity
                    //calculate (v n) n/A
                    Scalar tmp = fluxVars.KmvpNormal();
                    velocity = fluxVars.face().normal;
                    velocity *= tmp;
                    velocity /= scvfArea;
                    velocity *= (upwindWeight_ / up.viscosity() +
                                 (1 - upwindWeight_)/ dn.viscosity());

                    // add surface area for weighting purposes
                    boxSurface[vertIIdx][1] += scvfArea;
                    boxSurface[vertJIdx][1] += scvfArea;

                    velocityY[vertJIdx] += velocity[1];
                    velocityY[vertIIdx] += velocity[1];
                }
            }
#endif
        }
#ifdef VELOCITY_OUTPUT
        // normalize the velocities at the vertices
        // calculate the bounding box of the grid view
        VertexIterator vIt = this->gridView_().template begin<dim>();
        const VertexIterator vEndIt = this->gridView_().template end<dim>();
        for (; vIt!=vEndIt; ++vIt)
        {
            int i = this->problem_().vertexMapper().map(*vIt);

            //use vertiacl faces for vx and horizontal faces for vy calculation
            velocityX[i] /= boxSurface[i][0];
            if (dim >= 2)
            {
                velocityY[i] /= boxSurface[i][1];
            }
            if (dim == 3)
            {
                velocityZ[i] /= boxSurface[i][2];
            }
        }
#endif
        writer.attachVertexData(pressure, "P");
        writer.attachVertexData(delp, "delp");
#ifdef VELOCITY_OUTPUT // check if velocity output is demanded
        writer.attachVertexData(velocityX, "Vx");
        writer.attachVertexData(velocityY, "Vy");
        if (dim > 2)
            writer.attachVertexData(velocityZ, "Vz");
#endif
        writer.attachVertexData(moleFrac0, "x_H2O");
        writer.attachVertexData(moleFrac1, "x_N2");
        writer.attachVertexData(massFrac0, "X_H2O");
        writer.attachVertexData(massFrac1, "X_N2");
//        writer.attachVertexData(delFrac, "delFrac_TRAIL");
        writer.attachVertexData(rho, "rho");
        writer.attachVertexData(mu, "mu");
        writer.attachCellData(rank, "process rank");
    }

private:
    Scalar upwindWeight_;
};
}

#include "1p2cpropertydefaults.hh"

#endif
