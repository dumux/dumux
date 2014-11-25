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
 * \brief This file contains the data which is required to calculate the
 *        component fluxes over a face of a finite volume.
 *
 * This means concentration gradients, diffusion coefficients, mass fractions, etc.
 * at the integration point.
 */
#ifndef DUMUX_STOKESNC_FLUX_VARIABLES_HH
#define DUMUX_STOKESNC_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/freeflow/stokes/stokesfluxvariables.hh>

namespace Dumux
{

/*!
 * \ingroup BoxStokesncModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains data which is required to
 *        calculate the component fluxes over a face of a finite
 *        volume for the compositional n component Stokes model.
 *
 * This means concentration gradients, diffusion coefficient, mass fractions, etc.
 * at the integration point of a SCV or boundary face.
 */
template <class TypeTag>
class StokesncFluxVariables : public StokesFluxVariables<TypeTag>
{
    typedef StokesFluxVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    
    //dimensions
    enum {	dim = GridView::dimension };
    //phase indices
    enum {  phaseIdx = Indices::phaseIdx };
    //component indices
	enum {	phaseCompIdx = Indices::phaseCompIdx,
            transportCompIdx = Indices::transportCompIdx };
    //number of components
	enum {	numComponents = Indices::numComponents };
    
	typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
	typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

public:
	//Constrcutor calls ParentType function
	StokesncFluxVariables(const Problem &problem,
                          const Element &element,
                          const FVElementGeometry &fvGeometry,
                          const int fIdx,
                          const ElementVolumeVariables &elemVolVars,
                          const bool onBoundary = false)
		: ParentType(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary)
    {
        calculateValues_(problem, element, elemVolVars);
    }
	
	/*!
     * \brief Return the molar density \f$ \mathrm{[mol/m^3]} \f$ at the integration point.
     */
    const Scalar molarDensity() const
    { return molarDensity_; }
	
	/*!
     * \brief Return the mass fraction of a transported component at the integration point.
     */
    const Scalar massFraction(int compIdx) const
    { return massFraction_[compIdx]; }

    /*!
     * \brief Return the molar diffusion coefficient of a transported component at the integration point.
     */
    const Scalar diffusionCoeff(int compIdx) const
    { return diffusionCoeff_[compIdx]; }

    /*!
     * \brief Return the gradient of the mole fraction at the integration point.
     */
    const DimVector &moleFractionGrad(int compIdx) const
    { return moleFractionGrad_[compIdx]; }

	/*!
     * \brief Return the eddy diffusivity (if implemented).
     */
    const Scalar eddyDiffusivity() const
    { return 0; }

protected:
    void calculateValues_(const Problem &problem,
                          const Element &element,
                          const ElementVolumeVariables &elemVolVars)
    {

		// loop over all components
		for (int compIdx=0; compIdx<numComponents; compIdx++){
            massFraction_[compIdx] = Scalar(0.0);
            diffusionCoeff_[compIdx] = Scalar(0.0);
            moleFractionGrad_[compIdx] = Scalar(0.0);
                
			if (phaseCompIdx!=compIdx) //no transport equation parameters needed for the mass balance
			{
                molarDensity_ = Scalar(0.0);

				// calculate gradients and secondary variables at IPs
				for (int scvIdx = 0;
					 scvIdx < this->fvGeometry_.numScv;
					 scvIdx++) // loop over vertices of the element
				{
            
					molarDensity_ += elemVolVars[scvIdx].molarDensity()*
						this->face().shapeValue[scvIdx];
					massFraction_[compIdx] += elemVolVars[scvIdx].massFraction(compIdx) *
						this->face().shapeValue[scvIdx];
					diffusionCoeff_[compIdx] += elemVolVars[scvIdx].diffusionCoeff(compIdx) *
						this->face().shapeValue[scvIdx];

					// the gradient of the mole fraction at the IP
					for (int dimIdx=0; dimIdx<dim; ++dimIdx)
					{
						moleFractionGrad_[compIdx] +=
							this->face().grad[scvIdx][dimIdx] *
							elemVolVars[scvIdx].moleFraction(compIdx);
					}
				}
							
				Valgrind::CheckDefined(molarDensity_);
				Valgrind::CheckDefined(massFraction_[compIdx]);
				Valgrind::CheckDefined(diffusionCoeff_[compIdx]);
				Valgrind::CheckDefined(moleFractionGrad_[compIdx]);
			}
		}
    }
	
	Scalar molarDensity_;  
	Scalar massFraction_[numComponents];
    Scalar diffusionCoeff_[numComponents];
    DimVector moleFractionGrad_[numComponents];
};

} // end namespace

#endif
