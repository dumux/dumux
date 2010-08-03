// $Id: 2p2cnifluxvars.hh 3736 2010-06-15 09:52:10Z lauser $
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This file contains the data which is required to calculate
 *        all fluxes (mass of components and energy) over a face of a finite volume.
 *
 * This means pressure, concentration and temperature gradients, phase
 * densities at the integration point, etc.
 */
#ifndef DUMUX_2P2CNI_FLUX_DATA_HH
#define DUMUX_2P2CNI_FLUX_DATA_HH

#include <dumux/common/math.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCNIModel
 * \brief This template class contains the data which is required to
 *        calculate all fluxes (mass of components and energy) over a face of a finite
 *        volume for the non-isothermal two-phase, two-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class TwoPTwoCNIFluxVars : public TwoPTwoCFluxVars<TypeTag>
{
    typedef TwoPTwoCFluxVars<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SecondaryVars)) SecondaryVars;

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSecondaryVars)) ElementSecondaryVars;

    enum {
        dim           = GridView::dimension,
        dimWorld      = GridView::dimensionworld,

        numPhases     = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldVector<CoordScalar, dimWorld> Vector;

public:
    TwoPTwoCNIFluxVars(const Problem &problem,
                       const Element &element,
                       const FVElementGeometry &elemGeom,
                       int scvfIdx,
                       const ElementSecondaryVars &elemDat)
        : ParentType(problem, element, elemGeom, scvfIdx, elemDat)
    {
        // calculate temperature gradient using finite element
        // gradients
        Vector temperatureGrad(0);
        Vector tmp(0.0);
        for (int vertIdx = 0; vertIdx < elemGeom.numVertices; vertIdx++)
        {
            tmp = elemGeom.subContVolFace[scvfIdx].grad[vertIdx];
            tmp *= elemDat[vertIdx].temperature();
            temperatureGrad += tmp;
        }

        // The soil calculates the actual heat flux vector
        problem.spatialParameters().matrixHeatFlux(tmp,
                                                   *this,
                                                   elemDat,
                                                   temperatureGrad,
                                                   element,
                                                   elemGeom,
                                                   scvfIdx);
        // project the heat flux vector on the face's normal vector
        normalMatrixHeatFlux_ = tmp*elemGeom.subContVolFace[scvfIdx].normal;
    }

    /*!
     * \brief The total heat flux \f$[J/s]\f$ due to heat conduction
     *        of the rock matrix over the sub-control volume's face in
     *        direction of the face normal.
     */
    Scalar normalMatrixHeatFlux() const
    { return normalMatrixHeatFlux_; }

private:
    Scalar normalMatrixHeatFlux_;
};

} // end namepace

#endif
