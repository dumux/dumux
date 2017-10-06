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
#include "localjacobian.hh"
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
class StokesModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

public:
    //! \copydoc ImplicitModel::addOutputVtkFields
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > VelocityField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);
        ScalarField &pn = *writer.allocateManagedBuffer(numVertices);
        ScalarField &delP = *writer.allocateManagedBuffer(numVertices);
        ScalarField &rho = *writer.allocateManagedBuffer(numVertices);
        ScalarField &mu = *writer.allocateManagedBuffer(numVertices);
        VelocityField &velocity = *writer.template allocateManagedBuffer<Scalar, dim> (numVertices);

        unsigned numElements = this->gridView_().size(0);
        ScalarField &rank = *writer.allocateManagedBuffer(numElements);

        FVElementGeometry fvGeometry;
        VolumeVariables volVars;
        ElementBoundaryTypes elemBcTypes;

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
    }
};
}

#include "propertydefaults.hh"

#endif






















/*!
 * \file
 *
 * \brief Base class for all models which use the linear elasticity model.
 *        Adaption of the fully implicit scheme to the linear elasticity model.
 */
#ifndef DUMUX_ELASTIC_MODEL_HH
#define DUMUX_ELASTIC_MODEL_HH

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup ElasticBoxModel
 * \brief Adaption of the fully implicit scheme to the Navier-Stokes-equation.
 *
 * This model implements a linear elastic solid using Hooke's law as
 * stress-strain relation and a quasi-stationary momentum balance equation:
 \f[
  \boldsymbol{\sigma} = 2\,\my\,\boldsymbol{\epsilon} - \text{p}\, \boldsymbol{I}.
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
 (\textbf{u \dot grad})\boldsymbol{u} - \text{div} \boldsymbol{\sigma} - \varrho {\textbf{\textit{f}}} = 0 \;,
 \f]
 *
 * The equation is discretized using a vertex-centered finite volume (box)
 * scheme as spatial discretization.
 *
 */

template<class TypeTag >
class ElasticModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MechanicalLaw = typename GET_PROP_TYPE(TypeTag, MechanicalLaw);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SecondaryVariables = typename GET_PROP_TYPE(TypeTag, SecondaryVariables);

    static const int dim = GridView::dimension;
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
        // create the required scalar fields
        unsigned numVert = this->gridView_().size(dim);
        unsigned numElements = this->gridView_().size(0);

        auto& ux = *writer.allocateManagedBuffer(numVert);
        auto& uy = *writer.allocateManagedBuffer(numVert);
        auto& uz = *writer.allocateManagedBuffer(numVert);
        auto& uxExact = *writer.allocateManagedBuffer(numVert);
        auto& uyExact = *writer.allocateManagedBuffer(numVert);
        auto& uzExact = *writer.allocateManagedBuffer(numVert);
        auto& sigmax = *writer.template allocateManagedBuffer<Scalar, dim>(numElements);
        auto& sigmay = *writer.template allocateManagedBuffer<Scalar, dim>(numElements);
        auto& sigmaz = *writer.template allocateManagedBuffer<Scalar, dim>(numElements);
        auto& rank = *writer.allocateManagedBuffer(numElements);

        // initialize stress fields
        for (unsigned int i = 0; i < numElements; ++i)
        {
            sigmax[i] = 0;
            if (dim > 1)
                sigmay[i] = 0;
            if (dim > 2)
                sigmaz[i] = 0;
        }

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
            ElementSolution elemSol(numLocalDofs);
            for (unsigned int i = 0; i < numLocalDofs; ++i)
            {
                // only proceed for vertex dofs
                if (fe.localCoefficients().localKey(i).codim() != dim)
                    continue;

                const auto dofIdxGlobal = localIndexSet.index(i);

                const auto dofSol = sol[dofIdxGlobal];
                ux[dofIdxGlobal] = dofSol[Indices::u(0)];
                if (dim >= 2)
                    uy[dofIdxGlobal] = dofSol[Indices::u(1)];
                if (dim >= 3)
                    uz[dofIdxGlobal] = dofSol[Indices::u(2)];

                const auto vertexPos = element.template subEntity</*codim=*/dim>(fe.localCoefficients().localKey(i).subEntity()).geometry().center();
                const auto exactSolAtPos = this->problem_().exactSolution(vertexPos);
                uxExact[dofIdxGlobal] = exactSolAtPos[Indices::u(0)];
                if (dim >= 2)
                    uyExact[dofIdxGlobal] = exactSolAtPos[Indices::u(1)];
                if (dim >= 3)
                    uzExact[dofIdxGlobal] = exactSolAtPos[Indices::u(2)];

                elemSol[i] = std::move(dofSol);
            }

            // obtain element geometry
            auto eg = element.geometry();

            // evaluate shape function data and secondary variables at the cell center
            IpData ipData(eg, eg.local(eg.center()), fe.localBasis());
            SecondaryVariables secVars;
            secVars.update(elemSol, this->problem_(), element, ipData);

            // get the lame parameters
            const auto& lameParams = this->problem_().spatialParams().lameParams(element, secVars.priVars());

            // compute the stress tensor and add to container
            auto sigma = MechanicalLaw::stressTensor(element, ipData, secVars, elemSol, lameParams);
            sigmax[eIdx] += sigma[0];
            if (dim >= 2)
                sigmay[eIdx] += sigma[1];
            if (dim == 3)
                sigmaz[eIdx] += sigma[2];
        }

        writer.attachDofData(ux, "ux", true);
        if (dim >= 2)
            writer.attachDofData(uy, "uy", true);
        if (dim == 3)
            writer.attachDofData(uz, "uz", true);
        writer.attachDofData(uxExact, "ux_exact", true);
        if (dim >= 2)
            writer.attachDofData(uyExact, "uy_exact", true);
        if (dim == 3)
            writer.attachDofData(uzExact, "uz_exact", true);
        writer.attachCellData(sigmax, "stress X", dim);
        if (dim >= 2)
        writer.attachCellData(sigmay, "stress Y", dim);
        if (dim == 3)
        writer.attachCellData(sigmaz, "stress Z", dim);
    }
};
}
#include "propertydefaults.hh"
#endif
