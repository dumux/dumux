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
 * \brief This file contains the diffusion module for the vertex data
 *        of the fully coupled two-phase N-component model
 */
#ifndef DUMUX_MPNC_DIFFUSION_VOLUME_VARIABLES_HH
#define DUMUX_MPNC_DIFFUSION_VOLUME_VARIABLES_HH

namespace Dumux {

template<class TypeTag, bool enableDiffusion>
class MPNCVolumeVariablesDiffusion
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;

    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };
    enum { lPhaseIdx = FluidSystem::lPhaseIdx };
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };

public:
    MPNCVolumeVariablesDiffusion()
    {}

    template <class FluidState, class ParameterCache>
    void update(FluidState &fluidState,
                ParameterCache &paramCache,
                const VolumeVariables &volVars,
                const Problem &problem)
    {
        Valgrind::SetUndefined(*this);

        // diffusion coefficents in liquid
        diffCoeffL_[0] = 0.0;
        for (int compIdx = 1; compIdx < numComponents; ++compIdx) {
            diffCoeffL_[compIdx] =
                FluidSystem::binaryDiffusionCoefficient(fluidState,
                                                        paramCache,
                                                        lPhaseIdx,
                                                        0,
                                                        compIdx);
        }
        Valgrind::CheckDefined(diffCoeffL_);

        // diffusion coefficents in gas
        for (int compIIdx = 0; compIIdx < numComponents; ++compIIdx) {
            diffCoeffG_[compIIdx][compIIdx] = 0;
            for (int compJIdx = compIIdx + 1; compJIdx < numComponents; ++compJIdx) {
                diffCoeffG_[compIIdx][compJIdx] =
                        FluidSystem::binaryDiffusionCoefficient(fluidState,
                                                                paramCache,                                                        gPhaseIdx,
                                                        compIIdx,
                                                        compJIdx);

                // fill the symmetric part of the diffusion coefficent
                // matrix
                diffCoeffG_[compJIdx][compIIdx] = diffCoeffG_[compIIdx][compJIdx];
            }
        }
        Valgrind::CheckDefined(diffCoeffG_);
    };


    Scalar diffCoeff(int phaseIdx, int compIIdx, int compJIdx) const
    {
        if (phaseIdx == gPhaseIdx)
            // TODO: tensorial diffusion coefficients
            return diffCoeffG_[compIIdx][compJIdx];

        int i = std::min(compIIdx, compJIdx);
        int j = std::max(compIIdx, compJIdx);
        if (i != 0)
            return 0;
        return diffCoeffL_[j];
    };

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    {
        Valgrind::CheckDefined(diffCoeffL_);
        Valgrind::CheckDefined(diffCoeffG_);
    }


protected:
    // the diffusion coefficients for the porous medium for the
    // liquid phase
    Scalar diffCoeffL_[numComponents];

    // the diffusion coefficients for the porous medium for the
    // gas phase
    Scalar diffCoeffG_[numComponents][numComponents];
};

// dummy class for the case where diffusion is disabled
template<class TypeTag>
class MPNCVolumeVariablesDiffusion<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    MPNCVolumeVariablesDiffusion()
    {}

    template <class FluidState, class ParameterCache>
    void update(FluidState &fluidState,
                ParameterCache &paramCache,
                const VolumeVariables &volVars,
                const Problem &problem)
    { };

    Scalar diffCoeffL(int compIdx) const
    { return 0; };

    Scalar diffCoeffG(int compIIdx, int compJIdx) const
    { return 0; };

    /*!
     * \brief If running under valgrind this produces an error message
     *        if some of the object's attributes is undefined.
     */
    void checkDefined() const
    { }
};

};

#endif // DUMUX_MPNC_DIFFUSION_VOLUME_VARIABLES_HH
