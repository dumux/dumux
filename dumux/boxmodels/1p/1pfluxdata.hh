//$Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Onur Dogan                                   *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 *        the flux of the fluid over a face of a finite volume.
 */
#ifndef DUMUX_1P_FLUX_DATA_HH
#define DUMUX_1P_FLUX_DATA_HH

#include <dumux/auxiliary/math.hh>

namespace Dune
{

/*!
 * \ingroup OnePBoxModel
 * \brief This template class contains the data which is required to
 *        calculate the flux of the fluid over a face of a
 *        finite volume for the one-phase model.
 */
template <class TypeTag>
class OnePFluxData
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))    Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData)) VertexData;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef std::vector<VertexData>                      VertexDataArray;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume             SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace         SCVFace;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePIndices)) Indices;

public:
    OnePFluxData(const Problem &problem,
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
        viscosityAtIP = 0;
        pressureGrad = Scalar(0);

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

            // the pressure gradient
            tmp = feGrad;
            tmp *= elemDat[idx].pressure;
            pressureGrad += tmp;

            // fluid density
            densityAtIP +=
                elemDat[idx].density*face->shapeValue[idx];
            // fluid viscosity
            viscosityAtIP +=
                elemDat[idx].viscosity*face->shapeValue[idx];
        }

        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        tmp = problem.gravity();
        tmp *= densityAtIP;

        pressureGrad -= tmp;
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
        K.mv(pressureGrad, vDarcy);  // vDarcy = K * grad p
        vDarcyNormal = vDarcy*face->normal;

        // set the upstream and downstream vertices
        upstreamIdx = face->i;
        downstreamIdx = face->j;
        if (vDarcyNormal > 0)
            std::swap(upstreamIdx, downstreamIdx);
    }

public:
    const FVElementGeometry &fvElemGeom;
    const SCVFace *face;
    const SCV     *insideSCV;
    const SCV     *outsideSCV;

    // gradients
    GlobalPosition pressureGrad;

    // density of the fluid at the integration point
    Scalar densityAtIP;
    // viscosity of the fluid at the integration point
    Scalar viscosityAtIP;

    // darcy velocity in direction of the face normal
    Scalar vDarcyNormal;

    // local index of the upwind vertex
    int upstreamIdx;
    // local index of the downwind vertex
    int downstreamIdx;
};

} // end namepace

#endif
