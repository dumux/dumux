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
 * \brief Adaptation of the box scheme to the n-component Stokes model.
 */
#ifndef DUMUX_STOKESNC_MODEL_HH
#define DUMUX_STOKESNC_MODEL_HH

#include <dune/common/version.hh>

#include <dumux/freeflow/stokes/stokesmodel.hh>

#include "stokesnclocalresidual.hh"
#include "stokesncproperties.hh"

namespace Dumux {
/*!
 * \ingroup BoxStokesncModel
 * \brief Adaptation of the BOX scheme to the compositional Stokes model.
 *
 * This model implements an isothermal n-component Stokes flow of a fluid
 * solving a momentum balance, a mass balance and a conservation equation for each
 * component. When using mole fractions naturally the densities represent molar
 * densites
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
 * Component mass balance equations:
 * \f[
 \frac{\partial \left(\varrho_g X_g^\kappa\right)}{\partial t}
 + \boldsymbol{\nabla} \boldsymbol{\cdot} \left( \varrho_g {\boldsymbol{v}}_g X_g^\kappa
 - D^\kappa_g \varrho_g \frac{M^\kappa}{M_g} \boldsymbol{\nabla} X_g^\kappa \right)
 - q_g^\kappa = 0
 * \f]
 *
 * This is discretized using a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method in time.
 */
template<class TypeTag>
class StokesncModel : public StokesModel<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {  dim = GridView::dimension,
			transportCompIdx = Indices::transportCompIdx,
			phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx),
			useMoles = GET_PROP_VALUE(TypeTag, UseMoles),
			numComponents = Indices::numComponents 
	};

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

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
        ScalarField &pN = *writer.allocateManagedBuffer(numVertices);
        ScalarField &delP = *writer.allocateManagedBuffer(numVertices);
		ScalarField &T = *writer.allocateManagedBuffer(numVertices);

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

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator endit = this->gridView_().template end<0>();
        for (; elemIt != endit; ++elemIt)
        {
            int idx = this->elementMapper().map(*elemIt);
            rank[idx] = this->gridView_().comm().rank();

            fvGeometry.update(this->gridView_(), *elemIt);
            elemBcTypes.update(this->problem_(), *elemIt, fvGeometry);

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            int numLocalVerts = elemIt->subEntities(dim);
#else
            int numLocalVerts = elemIt->template count<dim>();
#endif
            for (int i = 0; i < numLocalVerts; ++i)
            {
                int vIdxGlobal = this->vertexMapper().map(*elemIt, i, dim);
                volVars.update(sol[vIdxGlobal],
                               this->problem_(),
                               *elemIt,
                               fvGeometry,
                               i,
                               false);

                pN[vIdxGlobal] = volVars.pressure()*scale_;
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
                velocity[vIdxGlobal] = volVars.velocity();
                velocity[vIdxGlobal] *= 1/scale_;
            }
        }
		writer.attachVertexData(T, "temperature");
        writer.attachVertexData(pN, "P");
        writer.attachVertexData(delP, "delP");
		
		if(useMoles) 
		{
			for (int j = 0; j < numComponents; ++j)
				{
					std::ostringstream oss;
					oss << "x"
					<< FluidSystem::componentName(j)
					<< FluidSystem::phaseName(phaseIdx);
					writer.attachVertexData(*moleFraction[j], oss.str().c_str());
				}
		} else {
			for (int j = 0; j < numComponents; ++j)
			{
				std::ostringstream oss;
				oss << "X"
				<< FluidSystem::componentName(j)
				<< FluidSystem::phaseName(phaseIdx);
				writer.attachVertexData(*massFraction[j], oss.str().c_str());
			}
		}
        
		
        writer.attachVertexData(rho, "rho");
        writer.attachVertexData(mu, "mu");
        writer.attachVertexData(velocity, "v", dim);
    }
};

}

#include "stokesncpropertydefaults.hh"

#endif
