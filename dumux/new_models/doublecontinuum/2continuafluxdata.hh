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
#ifndef DUMUX_2CONTINUA_FLUX_DATA_HH
#define DUMUX_2CONTINUA_FLUX_DATA_HH

#include <dumux/auxiliary/math.hh>

namespace Dune
{

/*!
 * \brief This template class contains the data which is required to
 *        calculate the fluxes of the fluid phases over a face of a
 *        finite volume for the double continuum, two-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class TwoContinuaFluxData
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
        numContinua = GET_PROP_VALUE(TypeTag, PTAG(NumContinua))
    };

    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume             SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace         SCVFace;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoContinuaIndices)) Indices;

public:
    TwoContinuaFluxData(const Problem &problem,
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

        for (int cont = 0; cont < numContinua; ++cont) {
            //densityAtIP[cont] = 0;
            viscosityAtIP[cont] = 0;
            pressureGrad[cont] = Scalar(0);
        }

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

            for (int cont = 0; cont < numContinua; ++cont) {
                // the pressure gradient
                tmp = feGrad;
                tmp *= elemDat[idx].pressure[cont];
                pressureGrad[cont] += tmp;

                // the concentration gradient
                tmp = feGrad;
                tmp *= elemDat[idx].molefraction[cont];
                concentrationGrad[cont] += tmp;

                // phase viscosity
                viscosityAtIP[cont] += elemDat[idx].viscosity[cont]*face->shapeValue[idx];
                }
        }

//        for (int cont = 0; cont < numContinua; ++cont) {
            // correct the pressure by the hydrostatic pressure due to
            // gravity
//            tmp = problem.gravity();
//            tmp *= densityAtIP[cont];
//           pressureGrad[cont] -= tmp;
//        }
    }

    void calculateVelocities_(const Problem &problem,
                              const Element &element,
                              const VertexDataArray &elemDat)
    {
        // calculate the permeability tensor. TODO: this should be
        // more flexible
        typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;

        Tensor K;
        GlobalPosition vDarcy;
        for (int cont = 0; cont < numContinua; ++cont) {
            const Tensor *Ki;
            const Tensor *Kj;
            if (cont == Indices::bloodCont) {
                Ki = &problem.soilBlood().K(insideSCV->global,
                                            element,
                                            insideSCV->local);
                Kj = &problem.soilBlood().K(outsideSCV->global,
                                            element,
                                            outsideSCV->local);
            }
            else if (cont == Indices::tissueCont) {
                Ki = &problem.soilTissue().K(insideSCV->global,
                                             element,
                                             insideSCV->local);
                Kj = &problem.soilTissue().K(outsideSCV->global,
                                             element,
                                             outsideSCV->local);
            }
            else
                continue;
            Dune::harmonicMeanMatrix(K, *Ki, *Kj);

            // temporary vector for the Darcy velocity
            K.mv(pressureGrad[cont], vDarcy);  // vDarcy = K * grad p
            vDarcyNormal[cont] = vDarcy*face->normal;

            // set the upstream and downstream vertices
            upstreamIdx[cont] = face->i;
            downstreamIdx[cont] = face->j;
            if (vDarcyNormal[cont] > 0)
                std::swap(upstreamIdx[cont], downstreamIdx[cont]);
        }
    }

    void calculateDiffCoeffPM_(const Problem &problem,
                               const Element &element,
                               const VertexDataArray &elemDat)
    {
        const VertexData &vDat_i = elemDat[face->i];
        const VertexData &vDat_j = elemDat[face->j];

        for (int cont = 0; cont < numContinua; ++cont) {
            // Diffusion coefficient in the porous medium
            diffCoeffPM[cont]
                = 1./2*(vDat_j.porosity[cont] * vDat_i.tortuosity[cont] * vDat_i.diffCoeff[cont] +
                        vDat_j.porosity[cont] * vDat_j.tortuosity[cont] * vDat_j.diffCoeff[cont]);
        }
    }

public:
    const FVElementGeometry &fvElemGeom;
    const SCVFace *face;
    const SCV     *insideSCV;
    const SCV     *outsideSCV;

    //! pressure gradients
    GlobalPosition pressureGrad[numContinua];
    //! concentration gradients
    GlobalPosition concentrationGrad[numContinua];

    //! darcy velocity in direction of the face normal
    Scalar vDarcyNormal[numContinua];

    //! densities of the fluid at the integration point
    //Scalar densityAtIP[numContinua];

    //! viscosity of the fluid at the integration point
    Scalar viscosityAtIP[numContinua];

    //! the effective diffusion coefficent in the porous medium
    Scalar diffCoeffPM[numContinua];

    //! local index of the upwind vertex for each phase
    int upstreamIdx[numContinua];
    //! local index of the downwind vertex for each phase
    int downstreamIdx[numContinua];
};

} // end namepace

#endif
