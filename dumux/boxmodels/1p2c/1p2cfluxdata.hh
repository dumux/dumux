// $Id: 1p2cfluxdata.hh 3838 2010-07-15 08:31:53Z bernd $
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

#include <dumux/common/math.hh>

namespace Dumux
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
    typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;

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
        molarDensityAtIP = 0;
        viscosityAtIP = 0;
        pressureGrad = Scalar(0);
        concentrationGrad = Scalar(0);

        calculateGradients_(problem, element, elemDat);
        calculateVelocities_(problem, element, elemDat);
        calculateDiffCoeffPM_(problem, element, elemDat);
        calculateDispersionTensor_(problem, element, elemDat);
    };

private:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const VertexDataArray &elemDat)
    {
        GlobalPosition tmp;
        if (!problem.spatialParameters().useTwoPointGradient(element, face->i, face->j)) {
            // use finite-element gradients
            tmp = 0.0;
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

                // molar phase density
                molarDensityAtIP +=
                    elemDat[idx].molarDensity*face->shapeValue[idx];

                // phase viscosity
                viscosityAtIP +=
                    elemDat[idx].viscosity*face->shapeValue[idx];
            }
        }
        else {
            // use two-point gradients
            tmp = element.geometry().corner(face->i);
            tmp -= element.geometry().corner(face->j);
            Scalar dist = tmp.two_norm();

            tmp  = face->normal;
            tmp /= face->normal.two_norm()*dist;

            pressureGrad       = tmp;
            pressureGrad      *= elemDat[face->j].pressure - elemDat[face->i].pressure;
            concentrationGrad  = tmp;
            concentrationGrad *= elemDat[face->j].molefraction - elemDat[face->i].molefraction;
            densityAtIP        = (elemDat[face->j].density + elemDat[face->i].density)/2;
            molarDensityAtIP        = (elemDat[face->j].molarDensity + elemDat[face->i].molarDensity)/2;
            viscosityAtIP      = (elemDat[face->j].viscosity  + elemDat[face->i].viscosity)/2;
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
        Tensor K;
        problem.spatialParameters().meanK(K,
                            problem.spatialParameters().intrinsicPermeability(element,
                                                                fvElemGeom,
                                                                face->i),
                            problem.spatialParameters().intrinsicPermeability(element,
                                                                fvElemGeom,
                                                                face->j));

        // vector for the Darcy velocity
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
            = 1./2*(vDat_i.porosity * vDat_i.tortuosity * vDat_i.diffCoeff +
                    vDat_j.porosity * vDat_j.tortuosity * vDat_j.diffCoeff);
    }

    void calculateDispersionTensor_(const Problem &problem,
            const Element &element,
            const VertexDataArray &elemDat)
    {
        const VertexData &vDat_i = elemDat[face->i];
        const VertexData &vDat_j = elemDat[face->j];

        //calculate dispersivity at the interface: [0]: alphaL = longitudinal disp. [m], [1] alphaT = transverse disp. [m]
        Dune::FieldVector<Scalar, 2> dispersivity(0);
        dispersivity[0] = 0.5 * (vDat_i.dispersivity[0] +  vDat_j.dispersivity[0]);
        dispersivity[1] = 0.5 * (vDat_i.dispersivity[1] +  vDat_j.dispersivity[1]);

        //calculate velocity at interface: v = -1/mu * vDarcy = -1/mu * K * grad(p)
        GlobalPosition velocity(vDarcy);
        velocity *= -1/viscosityAtIP;

        dispersionTensor = 0;

        //matrix multiplication of the velocity at the interface: vv^T
        for (int i=0; i<dim; i++)
        {
            for (int j = 0; j<dim; j++)
            {
                dispersionTensor[i][j]=velocity[i]*velocity[j];
            }
        }
        //normalize velocity product --> vv^T/||v||, [m/s]
        Scalar vNorm = velocity.two_norm();

        dispersionTensor /= vNorm;
        if (vNorm == 0)
        {
            dispersionTensor = 0;
        }

        //multiply with dispersivity difference: vv^T/||v||*(alphaL - alphaT), [m^2/s] --> alphaL = longitudinal disp., alphaT = transverse disp.
        dispersionTensor *= (dispersivity[0] - dispersivity[1]);

        //add ||v||*alphaT to the main diagonal:vv^T/||v||*(alphaL - alphaT) + ||v||*alphaT, [m^2/s]
        for (int i = 0; i<dim; i++)
        {
            dispersionTensor[i][i]+=vNorm*dispersivity[1];
        }
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

    //! molar density of the fluid at the integration point
    Scalar molarDensityAtIP;

    //! viscosity of the fluid at the integration point
    Scalar viscosityAtIP;

    //! the effective diffusion coefficent in the porous medium
    Scalar diffCoeffPM;

    //! the dispersion tensor in the porous medium
    Tensor dispersionTensor;

    //! darcy velocity at the face
   GlobalPosition vDarcy;

    //! darcy velocity in direction of the face normal
    Scalar vDarcyNormal;

    //! local index of the upwind vertex for each phase
    int upstreamIdx;
    //! local index of the downwind vertex for each phase
    int downstreamIdx;
};

} // end namepace

#endif
