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
#ifndef DUMUX_STOKES2C_FLUX_VARIABLES_HH
#define DUMUX_STOKES2C_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/freeflow/stokes/stokesfluxvariables.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Stokes2cIndices); //!< Enumerations for the compositional stokes models
}

/*!
 * \ingroup BoxStokes2cModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains data which is required to
 *        calculate the component fluxes over a face of a finite
 *        volume for the compositional Stokes model.
 *
 * This means concentration gradients, diffusion coefficient, mass fractions, etc.
 * at the integration point of a SCV or boundary face.
 */
template <class TypeTag>
class Stokes2cFluxVariables : public StokesFluxVariables<TypeTag>
{
    typedef StokesFluxVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension };
    enum { phaseIdx = Indices::phaseIdx };
    enum { transportCompIdx = Indices::transportCompIdx };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    Stokes2cFluxVariables(const Problem &problem,
                          const Element &element,
                          const FVElementGeometry &fvGeometry,
                          const int faceIdx,
                          const ElementVolumeVariables &elemVolVars,
                          const bool onBoundary = false)
        : ParentType(problem, element, fvGeometry, faceIdx, elemVolVars, onBoundary)
    {
        calculateValues_(problem, element, elemVolVars);
    }

    /*!
     * \brief Return the mass fraction of the transported component at the integration point.
     */
    const Scalar massFraction() const
    { return massFraction_; }

    /*!
     * \brief Return the mass fraction at the integration point.
     */
    DUNE_DEPRECATED_MSG("use massFraction() instead")
    Scalar massFractionAtIP() const
    { return massFraction(); }

    /*!
     * \brief Return the molar diffusion coefficient at the integration point.
     */
    const Scalar diffusionCoeff() const
    { return diffusionCoeff_; }

    /*!
     * \brief Return the eddy diffusivity (if implemented).
     */
    const Scalar eddyDiffusivity() const
    { return 0; }

    /*!
     * \brief Return the molar diffusion coefficient at the integration point.
     */
    DUNE_DEPRECATED_MSG("use diffusionCoeff() instead")
    Scalar diffusionCoeffAtIP() const
    { return diffusionCoeff(); }

    /*!
     * \brief Return the gradient of the mole fraction at the integration point.
     */
    const DimVector &moleFractionGrad() const
    { return moleFractionGrad_; }

    /*!
     * \brief Return the gradient of the mole fraction at the integration point.
     */
    DUNE_DEPRECATED_MSG("use moleFractionGrad() instead")
    const DimVector &moleFractionGradAtIP() const
    { return moleFractionGrad(); }


protected:
    void calculateValues_(const Problem &problem,
                          const Element &element,
                          const ElementVolumeVariables &elemVolVars)
    {
        massFraction_ = Scalar(0);
        diffusionCoeff_ = Scalar(0);
        moleFractionGrad_ = Scalar(0);

        // calculate gradients and secondary variables at IPs
        for (int idx = 0;
             idx < this->fvGeometry_.numVertices;
             idx++) // loop over vertices of the element
        {
            massFraction_ += elemVolVars[idx].fluidState().massFraction(phaseIdx, transportCompIdx) *
                this->face().shapeValue[idx];
            diffusionCoeff_ += elemVolVars[idx].diffusionCoeff() *
                this->face().shapeValue[idx];

            // the gradient of the mass fraction at the IP
            for (int dimIdx=0; dimIdx<dim; ++dimIdx)
            {
                moleFractionGrad_ +=
                    this->face().grad[idx][dimIdx] *
                    elemVolVars[idx].fluidState().moleFraction(phaseIdx, transportCompIdx);
            }
        }

        Valgrind::CheckDefined(massFraction_);
        Valgrind::CheckDefined(diffusionCoeff_);
        Valgrind::CheckDefined(moleFractionGrad_);
    }

    Scalar massFraction_;
    Scalar diffusionCoeff_;
    DimVector moleFractionGrad_;
};

} // end namespace

#endif
