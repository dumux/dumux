// $Id:$
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 *        all fluxes of fluid phases over a face of a finite volume.
 *
 * This means pressure and temperature gradients, phase densities at
 * the integration point, etc.
 */
#ifndef DUMUX_1P2C_FLUX_DATA_HH
#define DUMUX_1P2C_FLUX_DATA_HH

#include <dumux/auxiliary/math.hh>

namespace Dune
{

/*!
 * \brief This template class contains the data which is required to
 *        calculate the fluxes of the fluid phases over a face of a
 *        finite volume for the two-phase model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class OnePTwoCFluxData
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
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };

    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume             SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace         SCVFace;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;

public:
    OnePTwoCFluxData(const Problem &problem,
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
        concentrationGrad = Scalar(0);

        calculateGradients_(problem, element, elemDat);
        calculateVelocities_(problem, element, elemDat);
        calculateDiffCoeffPM_(problem, element, elemDat);
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

            // the concentration gradient
            tmp = feGrad;
            tmp *= elemDat[idx].molefraction;
            concentrationGrad += tmp;

            // phase density
            densityAtIP +=
                elemDat[idx].density*face->shapeValue[idx];

            // phase viscosity
            viscosityAtIP +=
                elemDat[idx].viscosity*face->shapeValue[idx];
        }

        // correct the pressure by the hydrostatic pressure due to
        // gravity
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity))) {
            tmp = problem.gravity();
            tmp *= densityAtIP;
            pressureGrad -= tmp;
        }
    }

    void calculateVelocities_(const Problem &problem,
                              const Element &element,
                              const VertexDataArray &elemDat)
    {
        // calculate the permeability tensor. TODO: this should be
        // more flexible
        typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;
        const Tensor &Ki = problem.soil().K(insideSCV->global,
                                            element,
                                            insideSCV->local);
        const Tensor &Kj = problem.soil().K(outsideSCV->global,
                                            element,
                                            outsideSCV->local);
        Tensor K;
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

    void calculateDiffCoeffPM_(const Problem &problem,
                               const Element &element,
                               const VertexDataArray &elemDat)
    {
        const VertexData &vDat_i = elemDat[face->i];
        const VertexData &vDat_j = elemDat[face->j];

        // Diffusion coefficient in the porous medium
        diffCoeffPM
            = 1./2*(vDat_i.density * vDat_i.porosity * vDat_i.tortuosity * vDat_i.diffCoeff +
                    vDat_j.density * vDat_j.porosity * vDat_j.tortuosity * vDat_j.diffCoeff);
    }

public:
    const FVElementGeometry &fvElemGeom;
    const SCVFace *face;
    const SCV     *insideSCV;
    const SCV     *outsideSCV;

    //! pressure gradient
    GlobalPosition pressureGrad;
    //! concentratrion gradient
    GlobalPosition concentrationGrad;

    //! density of the fluid at the integration point
    Scalar densityAtIP;

    //! viscosity of the fluid at the integration point
    Scalar viscosityAtIP;

    //! the effective diffusion coefficent in the porous medium
    Scalar diffCoeffPM;

    //! darcy velocity in direction of the face normal
    Scalar vDarcyNormal;

    //! local index of the upwind vertex for each phase
    int upstreamIdx;
    //! local index of the downwind vertex for each phase
    int downstreamIdx;
};

} // end namepace

#endif
