//$Id$
/*****************************************************************************
 *   Copyright (C) 2008,2009 by Klaus Mosthaf,                               *
 *                              Andreas Lauser,                              *
 *                              Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: klaus.mosthaf _at_ iws.uni-stuttgart.de                          *
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
 * This means pressure and temperature gradients, phase
 * densities at the intergration point, etc.
 */
#ifndef DUMUX_2P_FLUX_DATA_HH
#define DUMUX_2P_FLUX_DATA_HH

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
template <class TwoPTraits,
          class ProblemT,
          class VertexData>
class TwoPFluxData
{
    typedef typename ProblemT::DomainTraits::Scalar  Scalar;
    typedef typename ProblemT::DomainTraits::Grid    Grid;
    typedef ProblemT                                 Problem;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef std::vector<VertexData>                  VertexDataArray;

    static const int dim = Grid::dimensionworld;
    static const int localDim = Grid::dimension;
    typedef Dune::FieldVector<Scalar, dim>                   GlobalPosition;
    typedef Dune::FieldVector<Scalar, localDim>              LocalPosition;

    typedef typename Dune::FVElementGeometry<Grid>           FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume     SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    static const int numPhases = TwoPTraits::numPhases;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

    typedef TwoPTraits Tr;

public:
    TwoPFluxData(const Problem &problem,
                     const Element &element,
                     const FVElementGeometry &elemGeom,
                     int faceIdx,
                     const VertexDataArray &elemDat)
        : fvElemGeom(elemGeom)
    {
        face = &fvElemGeom.subContVolFace[faceIdx];

        int i = face->i;
        int j = face->j;

        insideSCV = &fvElemGeom.subContVol[i];
        outsideSCV = &fvElemGeom.subContVol[j];

        densityAtIP = 0;
        for (int phase = 0; phase < numPhases; ++phase) {
            pressureGrad[phase] = Scalar(0);
        }

        calculateGradients_(problem, element, elemDat);
        calculateVelocities_(problem, element, elemDat);
    };

private:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const VertexDataArray &elemDat)
    {
        // calculate gradients
        GlobalPosition tmp(0.0);
        for (int idx = 0;
             idx < fvElemGeom.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const LocalPosition &feGrad = face->grad[idx];

            // compute sum of pressure gradients for each phase
            for (int phase = 0; phase < numPhases; phase++)
            {
                // the pressure gradient
                tmp = feGrad;
                tmp *= elemDat[idx].pressure[phase];
                pressureGrad[phase] += tmp;

                // phase density
                densityAtIP[phase]
                    +=
                    elemDat[idx].density[phase] *
                    face->shapeValue[idx];
            }
        }

        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        for (int phase=0; phase < numPhases; phase++)
        {
            tmp = problem.gravity();
            tmp *= densityAtIP[phase];

            pressureGrad[phase] -= tmp;
        }
    }

    void calculateVelocities_(const Problem &problem,
                              const Element &element,
                              const VertexDataArray &elemDat)
    {
        // calculate the permeability tensor. TODO: this should be
        // more flexible
        typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;
        Tensor K;
        const Tensor &Ki = problem.soil().K(insideSCV->global,
                                            element,
                                            insideSCV->local);
        const Tensor &Kj = problem.soil().K(outsideSCV->global,
                                            element,
                                            outsideSCV->local);
        Dune::harmonicMeanMatrix(K, Ki, Kj);

        // temporary vector for the Darcy velocity
        GlobalPosition vDarcy;
        for (int phase=0; phase < numPhases; phase++)
        {
            K.mv(pressureGrad[phase], vDarcy);  // vDarcy = K * grad p
            vDarcyNormal[phase] = vDarcy*face->normal;
        }

        // set the upstream and downstream vertices
        for (int phase = 0; phase < numPhases; ++phase)
        {
            upstreamIdx[phase] = face->i;
            downstreamIdx[phase] = face->j;

            if (vDarcyNormal[phase] > 0) {
                std::swap(upstreamIdx[phase],
                          downstreamIdx[phase]);
            }
        }
    }

public:
    const FVElementGeometry &fvElemGeom;
    const SCVFace *face;
    const SCV     *insideSCV;
    const SCV     *outsideSCV;

    // gradients
    GlobalPosition pressureGrad[numPhases];

    // density of each face at the integration point
    PhasesVector densityAtIP;

    // darcy velocity in direction of the face normal
    PhasesVector vDarcyNormal;

    // local index of the upwind vertex for each phase
    int upstreamIdx[numPhases];
    // local index of the downwind vertex for each phase
    int downstreamIdx[numPhases];
};

} // end namepace

#endif
