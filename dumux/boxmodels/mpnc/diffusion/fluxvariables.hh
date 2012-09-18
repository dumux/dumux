// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
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

#include "../mpncproperties.hh"

namespace Dumux {

template<class TypeTag, bool enableDiffusion>
class MPNCFluxVariablesDiffusion
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    enum{dim = GridView::dimension};
    enum{numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum{numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum{wPhaseIdx = FluidSystem::wPhaseIdx};
    enum{nPhaseIdx = FluidSystem::nPhaseIdx};

    typedef Dune::FieldVector<Scalar, dim>  DimVector;

public:
    MPNCFluxVariablesDiffusion()
    {}

    void update(const Problem & problem,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const SCVFace & face,
                const ElementVolumeVariables & elemVolVars)
    {
        const unsigned int i = face.i;
        const unsigned int j = face.j;


        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx){
            for (int compIdx = 0; compIdx < numComponents; ++compIdx){
                moleFraction_[phaseIdx][compIdx] = 0. ;
                moleFractionGrad_[phaseIdx][compIdx] = 0. ;
            }
        }


        DimVector tmp ;
        for (int idx = 0;
             idx < fvGeometry.numFAP;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const DimVector & feGrad = face.grad[idx];

            // index for the element volume variables
            int volVarsIdx = face.fapIndices[idx];

            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++){
                for (int compIdx = 0; compIdx < numComponents; ++compIdx){

                    // calculate mole fractions at the integration points of the face
                    moleFraction_[phaseIdx][compIdx] += elemVolVars[volVarsIdx].fluidState().moleFraction(phaseIdx, compIdx)*
                    face.shapeValue[idx];

                    // calculate mole fraction gradients
                    tmp =  feGrad;
                    tmp *= elemVolVars[volVarsIdx].fluidState().moleFraction(phaseIdx, compIdx);
                    moleFractionGrad_[phaseIdx][compIdx] += tmp;
                }
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
            if (elemVolVars[i].fluidState().saturation(phaseIdx) <= 1e-4 ||
                elemVolVars[j].fluidState().saturation(phaseIdx) <= 1e-4)
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
                elemVolVars[i].fluidState().saturation(phaseIdx)/elemVolVars[i].porosity() *
                pow(elemVolVars[i].porosity() * elemVolVars[i].fluidState().saturation(phaseIdx), 7.0/3);
            Scalar red_j =
                elemVolVars[j].fluidState().saturation(phaseIdx)/elemVolVars[j].porosity() *
                pow(elemVolVars[j].porosity() * elemVolVars[j].fluidState().saturation(phaseIdx), 7.0/3);

            if (phaseIdx == wPhaseIdx) {
                // Liquid phase diffusion coefficients in the porous medium
                for (int i = 0; i < numComponents; ++i) {
                    // -> arithmetic mean
                    porousDiffCoeffL_[i]
                        = 1./2*(red_i * elemVolVars[i].diffCoeff(wPhaseIdx, 0, i) +
                                red_j * elemVolVars[j].diffCoeff(wPhaseIdx, 0, i));
                }
            }
            else {
                // Gas phase diffusion coefficients in the porous medium
                for (int i = 0; i < numComponents; ++i) {
                    for (int j = 0; j < numComponents; ++j) {
                        // -> arithmetic mean
                        porousDiffCoeffG_[i][j]
                            = 1./2*(red_i * elemVolVars[i].diffCoeff(nPhaseIdx, i, j) +
                                    red_j * elemVolVars[j].diffCoeff(nPhaseIdx, i, j));
                    }
                }
            }
        }
    }

    Scalar porousDiffCoeffL(const unsigned int compIdx) const
    {
        // TODO: tensorial diffusion coefficients
        return porousDiffCoeffL_[compIdx];
    }

    Scalar porousDiffCoeffG(const unsigned int compIIdx,
                            const unsigned int compJIdx) const
    {
        // TODO: tensorial diffusion coefficients
        return porousDiffCoeffG_[compIIdx][compJIdx];
    }

    Scalar moleFraction(const unsigned int phaseIdx,
                        const unsigned int compIdx) const
    { return moleFraction_[phaseIdx][compIdx]; }

    DUNE_DEPRECATED_MSG("use moleFraction() instead")
    Scalar moleFrac(const unsigned int phaseIdx,
                    const unsigned int compIdx) const
    { return moleFraction(phaseIdx, compIdx);}

    const DimVector &moleFractionGrad(const unsigned int phaseIdx,
                                  const unsigned int compIdx) const
    { return moleFractionGrad_[phaseIdx][compIdx];}

    DUNE_DEPRECATED_MSG("use moleFractionGrad() instead")
    const DimVector &moleFracGrad(const unsigned int phaseIdx,
                                  const unsigned int compIdx) const
    { return moleFractionGrad(phaseIdx, compIdx);}

protected:
    // the diffusion coefficients for the porous medium for the
    // liquid phase
    Scalar porousDiffCoeffL_[numComponents];

    // the diffusion coefficients for the porous medium for the
    // gas phase
    Scalar porousDiffCoeffG_[numComponents][numComponents];

    // the concentration gradients of all components in all phases
    DimVector moleFractionGrad_[numPhases][numComponents];

    // the mole fractions of each component at the integration point
    Scalar moleFraction_[numPhases][numComponents];
};


template<class TypeTag>
class MPNCFluxVariablesDiffusion<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    MPNCFluxVariablesDiffusion()
    {}

    void update(const Problem & problem,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const SCVFace & face,
                const ElementVolumeVariables & elemVolVars)
    {
    }
};

}

#endif // DUMUX_MPNC_DIFFUSION_FLUX_VARIABLES_HH
