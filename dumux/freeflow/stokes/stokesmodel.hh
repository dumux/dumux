// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2010-2012 by Katherina Baber, Klaus Mosthaf               *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
#ifndef DUMUX_STOKES_MODEL_HH
#define DUMUX_STOKES_MODEL_HH

#include <dumux/boxmodels/common/boxmodel.hh>

#include "stokeslocalresidual.hh"
#include "stokesnewtoncontroller.hh"
#include "stokeslocaljacobian.hh"
#include "stokesproblem.hh"
#include "stokesproperties.hh"

namespace Dumux
{
/*!
 * \ingroup BoxStokesModel
 * \brief Adaption of the box scheme to the Stokes model.
 *
 * This model implements laminar Stokes flow of a single fluid, solving a momentum balance:
 * \f[
\frac{\partial \left(\varrho_g {\boldsymbol{v}}_g\right)}{\partial t}
+ \boldsymbol{\nabla} \boldsymbol{\cdot} \left(p_g {\bf {I}}
- \mu_g \left(\boldsymbol{\nabla} \boldsymbol{v}_g
+ \boldsymbol{\nabla} \boldsymbol{v}_g^T\right)\right)
- \varrho_g {\bf g} = 0,
 * \f]
 *
 * and the mass balance equation:
 * \f[
\frac{\partial \varrho_g}{\partial t} + \boldsymbol{\nabla}\boldsymbol{\cdot}\left(\varrho_g {\boldsymbol{v}}_g\right) - q_g = 0
 * \f]
 *
 * This is discretized by a fully-coupled vertex-centered finite volume
 * (box) scheme in space and by the implicit Euler method in time.
 */
template<class TypeTag>
class StokesModel : public BoxModel<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

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
    void calculateFluxAcrossLayer(const SolutionVector &globalSol, Dune::FieldVector<Scalar, numEq> &flux, int axis, Scalar coordVal)
    {
        GlobalPosition globalI, globalJ;
        PrimaryVariables tmpFlux(0.0);
        int sign;

        FVElementGeometry fvElemGeom;
        ElementVolumeVariables elemVolVars;

        // Loop over elements
        ElementIterator elemIt = this->problem_.gridView().template begin<0>();
        ElementIterator endit = this->problem_.gridView().template end<0>();
        for (; elemIt != endit; ++elemIt)
        {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue;

            fvElemGeom.update(this->gridView_(), *elemIt);
            elemVolVars.update(this->problem_(), *elemIt, fvElemGeom);
            this->localResidual().evalFluxes(*elemIt, elemVolVars);

            bool hasLeft = false;
            bool hasRight = false;
            for (int i = 0; i < fvElemGeom.numVertices; i++) {
                const GlobalPosition &global = fvElemGeom.subContVol[i].global;
                if (globalI[axis] < coordVal)
                    hasLeft = true;
                else if (globalI[axis] >= coordVal)
                    hasRight = true;
            }
            if (!hasLeft || !hasRight)
                continue;

            for (int i = 0; i < fvElemGeom.numVertices; i++) {
                const int &globalIdx =
                    this->vertexMapper().map(*elemIt, i, dim);
                const GlobalPosition &global = fvElemGeom.subContVol[i].global;
                if (globalI[axis] < coordVal)
                    flux += this->localResidual().residual(i);
            }
        }

        flux = this->problem_.gridView().comm().sum(flux);
    }

    //! \copydoc BoxModel::addOutputVtkFields
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > VelocityField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);
        ScalarField &pN = *writer.allocateManagedBuffer(numVertices);
        ScalarField &delP = *writer.allocateManagedBuffer(numVertices);
        ScalarField &rho = *writer.allocateManagedBuffer(numVertices);
        ScalarField &mu = *writer.allocateManagedBuffer(numVertices);
        VelocityField &velocity = *writer.template allocateManagedBuffer<Scalar, dim> (numVertices);

        unsigned numElements = this->gridView_().size(0);
        ScalarField &rank = *writer.allocateManagedBuffer(numElements);

        FVElementGeometry fvElemGeom;
        VolumeVariables volVars;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator endit = this->gridView_().template end<0>();
        for (; elemIt != endit; ++elemIt)
        {
            int idx = this->elementMapper().map(*elemIt);
            rank[idx] = this->gridView_().comm().rank();

            fvElemGeom.update(this->gridView_(), *elemIt);
            elemBcTypes.update(this->problem_(), *elemIt, fvElemGeom);

            int numLocalVerts = elemIt->template count<dim>();
            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->vertexMapper().map(*elemIt, i, dim);
                volVars.update(sol[globalIdx],
                               this->problem_(),
                               *elemIt,
                               fvElemGeom,
                               i,
                               false);

                pN[globalIdx] = volVars.pressure();
                delP[globalIdx] = volVars.pressure() - 1e5;
                rho[globalIdx] = volVars.density();
                mu[globalIdx] = volVars.viscosity();
                velocity[globalIdx] = volVars.velocity();
            };
        }
        writer.attachVertexData(pN, "P");
        writer.attachVertexData(delP, "delP");
        writer.attachVertexData(rho, "rho");
        writer.attachVertexData(mu, "mu");
        writer.attachVertexData(velocity, "v", dim);
    }
};
}

#include "stokespropertydefaults.hh"

#endif
