// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Katherina Baber, Klaus Mosthaf                    *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
 * This means concentration gradients, diffusion coefficients, mass fractions, etc.
 * at the integration point.
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
    typedef typename GET_PROP_TYPE(TypeTag, Stokes2cIndices) Indices;

    enum { dim = GridView::dimension };
    enum { lCompIdx = Indices::lCompIdx };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIndex) };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> ScalarGradient;


public:
    Stokes2cFluxVariables(const Problem &problem,
                          const Element &element,
                          const FVElementGeometry &elemGeom,
                          int faceIdx,
                          const ElementVolumeVariables &elemVolVars,
                          bool onBoundary = false)
        : ParentType(problem, element, elemGeom, faceIdx, elemVolVars, onBoundary)
    {
        calculateValues_(problem, element, elemVolVars);
    }

    Scalar massFractionAtIP() const
    { return massFractionAtIP_; }

    Scalar diffusionCoeffAtIP() const
    { return diffusionCoeffAtIP_; }

    const ScalarGradient &massFractionGradAtIP() const
    { return massFractionGradAtIP_; }

protected:
    void calculateValues_(const Problem &problem,
                          const Element &element,
                          const ElementVolumeVariables &elemVolVars)
    {
        massFractionAtIP_ = Scalar(0);
        diffusionCoeffAtIP_ = Scalar(0);
        massFractionGradAtIP_ = Scalar(0);

        // calculate gradients and secondary variables at IPs
        for (int idx = 0;
             idx < this->fvGeom_.numVertices;
             idx++) // loop over vertices of the element
        {
            massFractionAtIP_ += elemVolVars[idx].fluidState().massFraction(phaseIdx, lCompIdx) *
                this->face().shapeValue[idx];
            diffusionCoeffAtIP_ += elemVolVars[idx].diffusionCoeff() *
                this->face().shapeValue[idx];

            // the gradient of the mass fraction at the IP
            for (int dimIdx=0; dimIdx<dim; ++dimIdx)
                massFractionGradAtIP_ +=
                    this->face().grad[idx][dimIdx]*
                    elemVolVars[idx].fluidState().massFraction(phaseIdx, lCompIdx);
        };

        Valgrind::CheckDefined(massFractionAtIP_);
        Valgrind::CheckDefined(diffusionCoeffAtIP_);
        Valgrind::CheckDefined(massFractionGradAtIP_);
    }

    Scalar massFractionAtIP_;
    Scalar diffusionCoeffAtIP_;
    ScalarGradient massFractionGradAtIP_;
};

} // end namespace

#endif
