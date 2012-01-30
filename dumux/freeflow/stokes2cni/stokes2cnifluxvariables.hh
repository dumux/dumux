// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Klaus Mosthaf                                     *
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
 *        energy fluxes over a face of a finite volume.
 *
 * This means concentration temperature gradients and heat conductivity at
 * integration point.
 */
#ifndef DUMUX_STOKES2CNI_FLUX_VARIABLES_HH
#define DUMUX_STOKES2CNI_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/freeflow/stokes2c/stokes2cfluxvariables.hh>

namespace Dumux
{

/*!
 * \ingroup BoxStokes2cniModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains data which is required to
 *        calculate the energy fluxes over a face of a finite
 *        volume for the non-isothermal compositional stokes model.
 *
 * This means temperature gradients, coefficients etc.
 * at the integration point.
 */
template <class TypeTag>
class Stokes2cniFluxVariables : public Stokes2cFluxVariables<TypeTag>
{
    typedef Stokes2cFluxVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    enum { dim = GridView::dimension };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> ScalarGradient;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    Stokes2cniFluxVariables(const Problem &problem,
                            const Element &element,
                            const FVElementGeometry &elemGeom,
                            int faceIdx,
                            const ElementVolumeVariables &elemVolVars,
                            bool onBoundary = false)
        : ParentType(problem, element, elemGeom, faceIdx, elemVolVars, onBoundary)
    {
        calculateValues_(problem, element, elemVolVars);
    }

    Scalar heatConductivityAtIP() const
    { return heatConductivityAtIP_; }

    const ScalarGradient &temperatureGradAtIP() const
    { return temperatureGradAtIP_; }

protected:
    void calculateValues_(const Problem &problem,
                          const Element &element,
                          const ElementVolumeVariables &elemVolVars)
    {
        heatConductivityAtIP_ = Scalar(0);
        temperatureGradAtIP_ = Scalar(0);

        // calculate gradients and secondary variables at IPs
        ScalarGradient tmp(0.0);
        for (int idx = 0;
             idx < this->fvGeom_.numVertices;
             idx++) // loop over vertices of the element
        {
            heatConductivityAtIP_ += elemVolVars[idx].heatConductivity() *
                this->face().shapeValue[idx];

            // the gradient of the temperature at the IP
            for (int dimIdx=0; dimIdx<dim; ++dimIdx)
                temperatureGradAtIP_ +=
                    this->face().grad[idx][dimIdx]*
                    elemVolVars[idx].temperature();
        }
        Valgrind::CheckDefined(heatConductivityAtIP_);
        Valgrind::CheckDefined(temperatureGradAtIP_);
    }

    Scalar heatConductivityAtIP_;
    ScalarGradient temperatureGradAtIP_;
};

} // end namespace

#endif
