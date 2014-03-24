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
/*!
 * \file
 *
 * \brief Adaptation of the box scheme to the two-component Stokes model.
 */
#ifndef DUMUX_STOKES2C_MODEL_HH
#define DUMUX_STOKES2C_MODEL_HH

#include <dumux/freeflow/stokes/stokesmodel.hh>

#include "stokes2clocalresidual.hh"
#include "stokes2cproperties.hh"

namespace Dumux {
/*!
 * \ingroup BoxStokes2cModel
 * \brief Adaptation of the BOX scheme to the compositional Stokes model.
 *
 * This model implements an isothermal two-component Stokes flow of a fluid
 * solving a momentum balance, a mass balance and a conservation equation for one component.
 *
 * Momentum Balance:
 * \f[
\frac{\partial \left(\varrho_g {\boldsymbol{v}}_g\right)}{\partial t}
+ \boldsymbol{\nabla} \boldsymbol{\cdot} \left(p_g {\bf {I}}
- \mu_g \left(\boldsymbol{\nabla} \boldsymbol{v}_g
+ \boldsymbol{\nabla} \boldsymbol{v}_g^T\right)\right)
- \varrho_g {\bf g} = 0,
 * \f]
 *
 * Mass balance equation:
 * \f[
\frac{\partial \varrho_g}{\partial t} + \boldsymbol{\nabla}\boldsymbol{\cdot}\left(\varrho_g {\boldsymbol{v}}_g\right) - q_g = 0
 * \f]
 *
 * Component mass balance equation:
 * \f[
 \frac{\partial \left(\varrho_g X_g^\kappa\right)}{\partial t}
 + \boldsymbol{\nabla} \boldsymbol{\cdot} \left( \varrho_g {\boldsymbol{v}}_g X_g^\kappa
 - D^\kappa_g \varrho_g \frac{M^\kappa}{M_g} \boldsymbol{\nabla} x_g^\kappa \right)
 - q_g^\kappa = 0
 * \f]
 *
 * This is discretized using a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as temporal discretization.
 */
    
template<class TypeTag>
    class Stokes2cModel : public StokesModel<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension };
    enum { transportCompIdx = Indices::transportCompIdx };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx) };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

public:
    /* Use Stokes n-component Model freeflow/stokesnc.The 2c model gets replaced by the more
     * general nc model. Adaption to the new model is straight forward. Please not that several
     * functions in the fluxvariables now need an input argument compIdx. Associated variables
     * are of size number of components. The nc model uses mole fraction formulations of the
     * transport equations as default. Mass fractions can be employed via the poperty UseMoles.*/
    DUNE_DEPRECATED_MSG("Use stokesnc model")
    Stokes2cModel(){}
    
    //! \copydoc ImplicitModel::addOutputVtkFields
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > VelocityField;

        const Scalar scale_ = GET_PROP_VALUE(TypeTag, Scaling);

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);
        ScalarField &pn = *writer.allocateManagedBuffer(numVertices);
        ScalarField &delP = *writer.allocateManagedBuffer(numVertices);
        ScalarField &Xw = *writer.allocateManagedBuffer(numVertices);
        ScalarField &rho = *writer.allocateManagedBuffer(numVertices);
        ScalarField &mu = *writer.allocateManagedBuffer(numVertices);
        VelocityField &velocity = *writer.template allocateManagedBuffer<Scalar, dim> (numVertices);

        unsigned numElements = this->gridView_().size(0);
        ScalarField &rank = *writer.allocateManagedBuffer(numElements);

        FVElementGeometry fvGeometry;
        VolumeVariables volVars;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator eEndIt = this->gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt)
        {
            int eIdx = this->elementMapper().map(*eIt);
            rank[eIdx] = this->gridView_().comm().rank();

            fvGeometry.update(this->gridView_(), *eIt);
            elemBcTypes.update(this->problem_(), *eIt, fvGeometry);

            int numLocalVerts = eIt->template count<dim>();
            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->vertexMapper().map(*eIt, i, dim);
                volVars.update(sol[globalIdx],
                               this->problem_(),
                               *eIt,
                               fvGeometry,
                               i,
                               false);

                pn[globalIdx] = volVars.pressure()*scale_;
                delP[globalIdx] = volVars.pressure()*scale_ - 1e5;
                Xw[globalIdx] = volVars.massFraction(phaseIdx, transportCompIdx);
                rho[globalIdx] = volVars.density()*scale_*scale_*scale_;
                mu[globalIdx] = volVars.dynamicViscosity()*scale_;
                velocity[globalIdx] = volVars.velocity();
                velocity[globalIdx] *= 1/scale_;
            }
        }
        writer.attachVertexData(pn, "P");
        writer.attachVertexData(delP, "delP");
        std::ostringstream outputNameX;
        outputNameX << "X^" << FluidSystem::componentName(transportCompIdx);
        writer.attachVertexData(Xw, outputNameX.str());
        writer.attachVertexData(rho, "rho");
        writer.attachVertexData(mu, "mu");
        writer.attachVertexData(velocity, "v", dim);
    }
};
}

#include "stokes2cpropertydefaults.hh"

#endif
