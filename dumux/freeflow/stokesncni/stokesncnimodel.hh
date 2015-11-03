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
 * \brief Adaption of the box scheme to the non-isothermal
 *        n-component Stokes model.
 */
#ifndef DUMUX_STOKESNCNI_MODEL_HH
#define DUMUX_STOKESNCNI_MODEL_HH

#include <dumux/freeflow/stokesnc/stokesncmodel.hh>

#include "stokesncnilocalresidual.hh"
#include "stokesncniproperties.hh"

namespace Dumux {
/*!
 * \ingroup BoxStokesncniModel
 * \brief Adaption of the box scheme to the non-isothermal
 *        n-component Stokes model.
 *
 * This model implements a non-isothermal n-component Stokes flow of a fluid
 * solving a momentum balance, a mass balance, conservation equations for \f$n-1\f$ components,
 * and one balance equation for the energy.
 *
 * The momentum balance:
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
 * The component mass balance equations:
 * \f[
 *    \frac{\partial \left(\varrho_g X_g^\kappa\right)}{\partial t}
 *    + \text{div} \left( \varrho_g {\boldsymbol{v}}_g X_g^\kappa
 *    - D^\kappa_g \varrho_g \frac{M^\kappa}{M_g} \textbf{grad}\, x_g^\kappa \right)
 *    - q_g^\kappa = 0
 * \f]
 * Please note that, even though it is n-component model, the diffusive
 * fluxes are still calculated with binary diffusion.
 *
 * The energy balance equation:
 * \f[
 *    \frac{\partial (\varrho_g  u_g)}{\partial t}
 *    + \text{div} \left( \varrho_g h_g {\boldsymbol{v}}_g
 *    - \sum_\kappa \left[ h^\kappa_g D^\kappa_g \varrho_g \frac{M^\kappa}{M_g}
 *                         \textbf{grad}\, x^\kappa_g \right]
 *    - \lambda_g \textbf{grad}\, T \right) - q_T = 0
 * \f]
 *
 * This is discretized by a fully-coupled vertex-centered finite volume
 * (box) scheme in space and by the implicit Euler method in time.
 */
template<class TypeTag>
class StokesncniModel : public StokesncModel<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension };
    enum { transportCompIdx = Indices::transportCompIdx };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx) };
    enum { useMoles = GET_PROP_VALUE(TypeTag, UseMoles) };
    enum { numComponents = Indices::numComponents };

   typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

public:
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
        ScalarField &T = *writer.allocateManagedBuffer(numVertices);
        ScalarField &h = *writer.allocateManagedBuffer(numVertices);
        ScalarField *moleFraction[numComponents];
        for (int i = 0; i < numComponents; ++i)
        moleFraction[i] = writer.template allocateManagedBuffer<Scalar, 1>(numVertices);

        ScalarField *massFraction[numComponents];
        for (int i = 0; i < numComponents; ++i)
        massFraction[i] = writer.template allocateManagedBuffer<Scalar, 1>(numVertices);

        ScalarField &rho = *writer.allocateManagedBuffer(numVertices);
        ScalarField &mu = *writer.allocateManagedBuffer(numVertices);
        VelocityField &velocity = *writer.template allocateManagedBuffer<Scalar, dim> (numVertices);

        unsigned numElements = this->gridView_().size(0);
        ScalarField &rank = *writer.allocateManagedBuffer(numElements);

        FVElementGeometry fvGeometry;
        VolumeVariables volVars;
        ElementBoundaryTypes elemBcTypes;

        for (const auto& element : Dune::elements(this->gridView_()))
        {
            int idx = this->elementMapper().index(element);

            rank[idx] = this->gridView_().comm().rank();

            fvGeometry.update(this->gridView_(), element);
            elemBcTypes.update(this->problem_(), element, fvGeometry);

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

                pn[vIdxGlobal] = volVars.pressure()*scale_;
                delP[vIdxGlobal] = volVars.pressure()*scale_ - 1e5;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    (*moleFraction[compIdx])[vIdxGlobal]= volVars.moleFraction(compIdx);
                    (*massFraction[compIdx])[vIdxGlobal]= volVars.massFraction(compIdx);
                    Valgrind::CheckDefined((*moleFraction[compIdx])[vIdxGlobal]);
                    Valgrind::CheckDefined((*massFraction[compIdx])[vIdxGlobal]);
                }

                T   [vIdxGlobal] = volVars.temperature();

                rho[vIdxGlobal] = volVars.density()*scale_*scale_*scale_;
                mu[vIdxGlobal] = volVars.dynamicViscosity()*scale_;
                h[vIdxGlobal] = volVars.enthalpy();
                velocity[vIdxGlobal] = volVars.velocity();
                velocity[vIdxGlobal] *= 1/scale_;
            }
        }
        writer.attachVertexData(T, "temperature");
        writer.attachVertexData(pn, "pg");
        writer.attachVertexData(delP, "delP");

        for (int j = 0; j < numComponents; ++j)
        {
            std::ostringstream moleFrac, massFrac;
            moleFrac << "x_" << FluidSystem::phaseName(phaseIdx)
                     << "^" << FluidSystem::componentName(j);
            writer.attachVertexData(*moleFraction[j], moleFrac.str().c_str());

            massFrac << "X_" << FluidSystem::phaseName(phaseIdx)
                     << "^" << FluidSystem::componentName(j);
            writer.attachVertexData(*massFraction[j], massFrac.str().c_str());
        }

        writer.attachVertexData(h, "h");
        writer.attachVertexData(rho, "rho");
        writer.attachVertexData(mu, "mu");
        writer.attachVertexData(velocity, "v", dim);

    }
};

}

#include "stokesncnipropertydefaults.hh"

#endif
