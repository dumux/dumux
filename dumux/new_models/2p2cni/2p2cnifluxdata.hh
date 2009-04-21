/*****************************************************************************
 *   Copyright (C) 2008,2009 by Melanie Darcis                               *
 *                              Klaus Mosthaf,                               *
 *                              Andreas Lauser,                              *
 *                              Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: melanie.darcis _at_ iws.uni-stuttgart.de                         *
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
 *        all fluxes of components over a face of a finite volume.
 *
 * This means pressure, concentration and temperature gradients, phase
 * densities at the intergration point, etc.
 */
#ifndef DUMUX_2P2CNI_FLUX_DATA_HH
#define DUMUX_2P2CNI_FLUX_DATA_HH

#include <dumux/auxiliary/math.hh>

namespace Dune
{

/*!
 * \brief This template class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the two-phase, two-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TwoPTwoCNITraits, 
          class ProblemT,
          class VertexData>
class TwoPTwoCNIFluxData : public TwoPTwoCFluxData<TwoPTwoCNITraits, ProblemT, VertexData>
{
    typedef typename ProblemT::DomainTraits::Scalar  Scalar;
    typedef ProblemT                                 Problem;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef std::vector<VertexData>                  VertexDataArray;

    static const int dim = Grid::dimensionworld;
    typedef Dune::FieldVector<Scalar, dim>           GlobalPosition;

    typedef TwoPTwoCNITraits Tr;
    typedef TwoPTwoCFluxData<TwoPTwoCNITraits, ProblemT, VertexData> ParentT;

public:
    TwoPTwoCFluxData(const Problem &problem,
                     const Element &element,
                     const FVElementGeometry &elemGeom,
                     int faceIdx,
                     const VertexDataArray &elemDat)
        : ParentType(problem, element, elemGeom, faceIdx, elemDat)
    {
        tempGrad = 0;

        // Harmonic mean of the heat conducitivities of the
        // sub control volumes adjacent to the face
        heatCondAtIp = harmonicMean(eDat.vertex[i].heatCond,
                                    eDat.vertex[j].heatCond);
        
        // calculate temperature gradient
        GlobalPosition tmp(0.0);
        for (int idx = 0; 
             idx < this->fvElemGeom.numVertices;
             idx++) // loop over adjacent vertices
        {
            tmp = face->grad[idx];
            tmp *= elemDat_[vertIdx].temperature;
            tempGrad += tmp;
        }
    }
      
    //! temperature gradient
    GlobalPosition temperatureGrad;
    
    //! heat conductivity at integration point
    Scalar heatContAtIp;
};

} // end namepace

#endif
