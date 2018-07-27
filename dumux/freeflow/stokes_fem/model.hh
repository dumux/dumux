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
#ifndef DUMUX_STOKES_MODEL_HH
#define DUMUX_STOKES_MODEL_HH

/*!
 * \file
 * \brief Base class for all models which use the Stokes box model.
 */

#include <dumux/implicit/model.hh>

#include "localresidual.hh"
#include "newtoncontroller.hh"
//#include "localjacobian.hh"
#include "problem.hh"
#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup BoxStokesModel
 * \brief Adaption of the box scheme to the Stokes model.
 *
 * This model implements laminar Stokes flow of a single fluid, solving the momentum balance equation
 * \f[
 *    \frac{\partial \left(\varrho_g {\boldsymbol{v}}_g\right)}{\partial t}
 *    + \text{div} \left( p_g {\bf {I}}
 *    - \mu_g \left( \textbf{grad}\, \boldsymbol{v}_g
 *                   + \textbf{grad}\, \boldsymbol{v}_g^T \right) \right)
 *    - \varrho_g {\bf g} = 0
 * \f]
 * By setting the property <code>EnableNavierStokes</code> to <code>true</code> the Navier-Stokes
 * equation can be solved. In this case an additional term
 * \f[
 *    + \text{div} \left( \varrho_g \boldsymbol{v}_g \boldsymbol{v}_g \right)
 * \f]
 * is added to the momentum balance equation.
 *
 * The mass balance equation:
 * \f[
 *    \frac{\partial \varrho_g}{\partial t}
 *    + \text{div} \left(\varrho_g {\boldsymbol{v}}_g\right) - q_g = 0
 * \f]
 *
 * This is discretized by a fully-coupled vertex-centered finite volume
 * (box) scheme in space and by the implicit Euler method in time.
 */
template<class TypeTag>
class StokesFemModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
/*
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

//    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
//    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
//    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
*/

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SecondaryVariables = typename GET_PROP_TYPE(TypeTag, SecondaryVariables);

    static const int dim = GridView::dimension;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;


public:
    /*!
     * \brief Calculate the fluxes across a certain layer in the domain.
     * The layer is situated perpendicular to the coordinate axis "coord" and cuts
     * the axis at the value "coordVal".
     *
     * \param globalSol The global solution vector
     * \param flux A vector to store the flux
     * \param axis The dimension, perpendicular to which the layer is situated
     * \param coordVal The (Scalar) coordinate on the axis, at which the layer is situated
     */
    /*
    void calculateFluxAcrossLayer(const SolutionVector &globalSol, Dune::FieldVector<Scalar, numEq> &flux, int axis, Scalar coordVal)
    {
        GlobalPosition globalI, globalJ;
        PrimaryVariables tmpFlux(0.0);

        FVElementGeometry fvGeometry;
        ElementVolumeVariables elemVolVars;

        // Loop over elements
        for (const auto& element : elements(this->problem_.gridView(), Dune::Partitions::interior))
        {
            fvGeometry.update(this->gridView_(), element);
            elemVolVars.update(this->problem_(), element, fvGeometry);
            this->localResidual().evalFluxes(element, elemVolVars);

            bool hasLeft = false;
            bool hasRight = false;
            for (int i = 0; i < fvGeometry.numScv; i++) {
                const GlobalPosition &globalPos = fvGeometry.subContVol[i].global;
                if (globalI[axis] < coordVal)
                    hasLeft = true;
                else if (globalI[axis] >= coordVal)
                    hasRight = true;
            }
            if (!hasLeft || !hasRight)
                continue;

            for (int i = 0; i < fvGeometry.numScv; i++) {
                const GlobalPosition &globalPos = fvGeometry.subContVol[i].global;
                if (globalI[axis] < coordVal)
                    flux += this->localResidual().residual(i);
            }
        }

        flux = this->problem_.gridView().comm().sum(flux);
    }
    */

    //! \copydoc ImplicitModel::addOutputVtkFields
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        //create the required scalar fields
    unsigned numVert = this->gridView_().size(dim);
    unsigned numElements = this->gridView_().size(0);

    auto& vx = *writer.allocateManagedBuffer(numVert);
    auto& vy = *writer.allocateManagedBuffer(numVert);
    auto& vz = *writer.allocateManagedBuffer(numVert);

    //exact solution
    auto& uxExact = *writer.allocateManagedBuffer(numVert);
    auto& uyExact = *writer.allocateManagedBuffer(numVert);
    auto& uzExact = *writer.allocateManagedBuffer(numVert);
    auto& pExact  = *writer.allocateManagedBuffer(numVert);
    auto& velocityExact = *writer.template allocateManagedBuffer<Scalar, dim> (numVert);

    auto& velocity = *writer.template allocateManagedBuffer<Scalar, dim> (numVert);

    auto& p = *writer.allocateManagedBuffer(numVert);

    //      auto& p = *writer.template allocateManagedBuffer<Scalar, dim>(numElements);

/*
    ScalarField &pn = *writer.allocateManagedBuffer(numVertices);
    ScalarField &delP = *writer.allocateManagedBuffer(numVertices);
    ScalarField &rho = *writer.allocateManagedBuffer(numVertices);
    ScalarField &mu = *writer.allocateManagedBuffer(numVertices);
    VelocityField &velocity = *writer.template allocateManagedBuffer<Scalar, dim> (numVertices);
*/

    auto& rank = *writer.allocateManagedBuffer(numElements);


    auto localView = this->feBasis().localView();
    auto localIndexSet = this->feBasis().localIndexSet();



    for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
    {
        auto eIdx = this->problem_().model().elementMapper().index(element);

        // rank output
        rank[eIdx] = this->gridView_().comm().rank();

        // bind local restrictions to element
        localView.bind(element);
        localIndexSet.bind(localView);

        // obtain the finite element
        const auto& fe = localView.tree().finiteElement();

        // loop over dofs inside the element and store solution at vertices
        // TODO HOW TO INCLUDE SUBSAMPLING
        const auto numLocalDofs = fe.localBasis().size();

        // obtain element geometry
        auto eg = element.geometry();
        // evaluate shape function data and secondary variables at the cell center
        IpData ipData(eg, eg.local(eg.center()), fe.localBasis());
        SecondaryVariables secVars;

        for (unsigned int i = 0; i < numLocalDofs; ++i)
        {
            // only proceed for vertex dofs
            if (fe.localCoefficients().localKey(i).codim() != dim)
                continue;

            const auto dofIdxGlobal = localIndexSet.index(i);
            const auto& dofSol = sol[dofIdxGlobal];

            vx[dofIdxGlobal] = dofSol[Indices::v(0)];
            velocity[dofIdxGlobal][Indices::v(0)] = dofSol[Indices::v(0)];
            if (dim >= 2){
               vy[dofIdxGlobal] = dofSol[Indices::v(1)];
               velocity[dofIdxGlobal][Indices::v(1)] = dofSol[Indices::v(1)];
            }
            if (dim >= 3){
               vz[dofIdxGlobal] = dofSol[Indices::v(2)];
               velocity[dofIdxGlobal][Indices::v(2)] = dofSol[Indices::v(2)];
            }

            //analytical solution
            const auto vertexPos = element.template subEntity</*codim=*/dim>(fe.localCoefficients().localKey(i).subEntity()).geometry().center();
            const auto exactSolAtPos = this->problem_().analyticalSolution(vertexPos);
            uxExact[dofIdxGlobal] = exactSolAtPos[Indices::v(0)];
            velocityExact[dofIdxGlobal][Indices::v(0)] = exactSolAtPos[Indices::v(0)];
            if (dim >= 2){
                uyExact[dofIdxGlobal] = exactSolAtPos[Indices::v(1)];
                velocityExact[dofIdxGlobal][Indices::v(1)] = exactSolAtPos[Indices::v(1)];
            }
            if (dim >= 3){
                uzExact[dofIdxGlobal] = exactSolAtPos[Indices::v(2)];
                velocityExact[dofIdxGlobal][Indices::v(2)] = exactSolAtPos[Indices::v(2)];
            }

            pExact[dofIdxGlobal] = exactSolAtPos[Indices::pressureIdx];



            // compute the stress tensor and add to container
            p[dofIdxGlobal] = dofSol[dim];
        }
    }


    writer.attachDofData(vx, "vx", true);
    writer.attachDofData(uxExact, "vxExact", true);
    if (dim >= 2){
        writer.attachDofData(vy, "vy", true);
        writer.attachDofData(uyExact, "vyExact", true);
    }
    if (dim == 3){
        writer.attachDofData(vz, "vz", true);
        writer.attachDofData(uzExact, "vzExact", true);
    }

    writer.attachDofData(p, "p", true);
    writer.attachDofData(pExact, "pExact", true);

    writer.attachVertexData(velocity, "v", dim);
    writer.attachVertexData(velocityExact, "vExact", dim);


/*        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > VelocityField;

    // create the required scalar fields
    unsigned numVertices = this->gridView_().size(dim);

    //added to resemble elastic
    unsigned numElements = this->gridView_().size(0);

    ScalarField &pn = *writer.allocateManagedBuffer(numVertices);
    ScalarField &delP = *writer.allocateManagedBuffer(numVertices);
    ScalarField &rho = *writer.allocateManagedBuffer(numVertices);
    ScalarField &mu = *writer.allocateManagedBuffer(numVertices);
    VelocityField &velocity = *writer.template allocateManagedBuffer<Scalar, dim> (numVertices);

    ScalarField &rank = *writer.allocateManagedBuffer(numElements);

    //        FVElementGeometry fvGeometry;
    //        VolumeVariables volVars;
    //        ElementBoundaryTypes elemBcTypes;

    for (const auto& element : elements(this->gridView_()))
    {
    int eIdx = this->elementMapper().index(element);

    rank[eIdx] = this->gridView_().comm().rank();

    fvGeometry.update(this->gridView_(), element);
    elemBcTypes.update(this->problem_(), element);

    int numLocalVerts = element.subEntities(dim);

    for (int i = 0; i < numLocalVerts; ++i)
    {
    int vIdxGlobal = this->vertexMapper().subIndex(element, i, dim);

    volVars.update(sol[vIdxGlobal],
                   this->problem_(),
                   element,
                   fvGeometry,
                   i,
                   false);

    pn[vIdxGlobal] = volVars.pressure();
    delP[vIdxGlobal] = volVars.pressure() - 1e5;
    rho[vIdxGlobal] = volVars.density();
    mu[vIdxGlobal] = volVars.dynamicViscosity();
    velocity[vIdxGlobal] = volVars.velocity();
    }
    }
    writer.attachVertexData(pn, "P");
    writer.attachVertexData(delP, "delP");
    writer.attachVertexData(rho, "rho");
    writer.attachVertexData(mu, "mu");
    writer.attachVertexData(velocity, "v", dim);
*/
    }
};


}

#include "propertydefaults.hh"

#endif
