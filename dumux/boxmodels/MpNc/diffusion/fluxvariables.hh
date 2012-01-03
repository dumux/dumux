// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief This file contains the diffusion module for the flux data of
 *        the fully coupled two-phase N-component model
 */
#ifndef DUMUX_MPNC_DIFFUSION_FLUX_VARIABLES_HH
#define DUMUX_MPNC_DIFFUSION_FLUX_VARIABLES_HH

#include <dune/common/fvector.hh>

#include "../MpNcproperties.hh"

namespace Dumux {

template<class TypeTag, bool enableDiffusion>
class MPNCFluxVariablesDiffusion
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;

    typedef typename GridView::template Codim<0>::Entity Element;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),

        lPhaseIdx = FluidSystem::lPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx,
    };

    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  Vector;
    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldVector<Scalar, numPhases>      PhasesVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;

public:
    MPNCFluxVariablesDiffusion()
    {}

    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvfIdx,
                const ElementVolumeVariables &vDat)
    {
        int i = elemGeom.subContVolFace[scvfIdx].i;
        int j = elemGeom.subContVolFace[scvfIdx].j;

        for (int phase = 0; phase < numPhases; ++phase) {
            for (int comp = 0; comp < numComponents; ++comp) {
                moleFrac_[phase][comp]  = vDat[i].fluidState().moleFraction(phase, comp);
                moleFrac_[phase][comp] += vDat[j].fluidState().moleFraction(phase, comp);
                moleFrac_[phase][comp] /= 2;
            }
        }

        // update the concentration gradients using two-point
        // gradients
        const GlobalPosition &normal = elemGeom.subContVolFace[scvfIdx].normal;

        GlobalPosition tmp = element.geometry().corner(j);
        tmp -= element.geometry().corner(i);
        Scalar dist = tmp.two_norm()*normal.two_norm();
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            // concentration gradients
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                moleFracGrad_[phaseIdx][compIdx] = normal;
                moleFracGrad_[phaseIdx][compIdx]
                    *=
                    (vDat[j].fluidState().moleFraction(phaseIdx, compIdx) -
                     vDat[i].fluidState().moleFraction(phaseIdx, compIdx))
                    / dist;
            }
        }

        // initialize the diffusion coefficients to zero
        for (int i = 0; i < numComponents; ++i) {
            porousDiffCoeffL_[i] = 0.0;
            for (int j = 0; j < numComponents; ++j)
                porousDiffCoeffG_[i][j] = 0.0;
        }

        // calculate the diffusion coefficients at the integration
        // point in the porous medium
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // make sure to only calculate diffusion coefficents
            // for phases which exist in both finite volumes
            if (vDat[i].fluidState().saturation(phaseIdx) <= 1e-4 ||
                vDat[j].fluidState().saturation(phaseIdx) <= 1e-4)
            {
                continue;
            }

            // reduction factor for the diffusion coefficients in the
            // porous medium in nodes i and j. this is the tortuosity
            // times porosity times phase saturation at the nodes i
            // and j
            //
            // TODO (?): move this calculation to the soil (possibly
            // that's a bad idea, though)
            Scalar red_i =
                vDat[i].fluidState().saturation(phaseIdx)/vDat[i].porosity() *
                pow(vDat[i].porosity() * vDat[i].fluidState().saturation(phaseIdx), 7.0/3);
            Scalar red_j =
                vDat[j].fluidState().saturation(phaseIdx)/vDat[j].porosity() *
                pow(vDat[j].porosity() * vDat[j].fluidState().saturation(phaseIdx), 7.0/3);

            if (phaseIdx == FluidSystem::lPhaseIdx) {
                // Liquid phase diffusion coefficients in the porous medium
                for (int i = 0; i < numComponents; ++i) {
                    // -> arithmetic mean
                    porousDiffCoeffL_[i]
                        = 1./2*(red_i * vDat[i].diffCoeff(lPhaseIdx, 0, i) +
                                red_j * vDat[j].diffCoeff(lPhaseIdx, 0, i));
                }
            }
            else {
                // Gas phase diffusion coefficients in the porous medium
                for (int i = 0; i < numComponents; ++i) {
                    for (int j = 0; j < numComponents; ++j) {
                        // -> arithmetic mean
                        porousDiffCoeffG_[i][j]
                            = 1./2*(red_i * vDat[i].diffCoeff(gPhaseIdx, i, j) +
                                    red_j * vDat[j].diffCoeff(gPhaseIdx, i, j));
                    }
                }
            }
        }
    };

    Scalar porousDiffCoeffL(int compIdx) const
    {
        // TODO: tensorial diffusion coefficients
        return porousDiffCoeffL_[compIdx];
    };

    Scalar porousDiffCoeffG(int compIIdx, int compJIdx) const
    {
        // TODO: tensorial diffusion coefficients
        return porousDiffCoeffG_[compIIdx][compJIdx];
    };

    Scalar moleFrac(int phaseIdx,
                    int compIdx) const
    {
        return moleFrac_[phaseIdx][compIdx];
    };

    const GlobalPosition &moleFracGrad(int phaseIdx,
                                            int compIdx) const
    {
        return moleFracGrad_[phaseIdx][compIdx];
    };

protected:
    // the diffusion coefficients for the porous medium for the
    // liquid phase
    Scalar porousDiffCoeffL_[numComponents];

    // the diffusion coefficients for the porous medium for the
    // gas phase
    Scalar porousDiffCoeffG_[numComponents][numComponents];

    // the concentration gradients of all components in all phases
    GlobalPosition moleFracGrad_[numPhases][numComponents];

    // the mole fractions of each component at the integration point
    Scalar moleFrac_[numPhases][numComponents];
};


template<class TypeTag>
class MPNCFluxVariablesDiffusion<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

public:
    MPNCFluxVariablesDiffusion()
    {}

    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvfIdx,
                const ElementVolumeVariables &vDat)
    {
    };
};

}

#endif // DUMUX_MPNC_DIFFUSION_FLUX_VARIABLES_HH
