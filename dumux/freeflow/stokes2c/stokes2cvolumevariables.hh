// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Katherina Baber, Klaus Mosthaf                    *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
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
 * \brief Contains the quantities which are constant within a
 *        finite volume in the compositional Stokes model.
 */
#ifndef DUMUX_STOKES2C_VOLUME_VARIABLES_HH
#define DUMUX_STOKES2C_VOLUME_VARIABLES_HH

#include <dumux/freeflow/stokes/stokesvolumevariables.hh>
#include "stokes2cproperties.hh"

namespace Dumux
{

/*!
 * \ingroup BoxStokes2cModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-component Stokes box model.
 */
template <class TypeTag>
class Stokes2cVolumeVariables : public StokesVolumeVariables<TypeTag>
{
    typedef StokesVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Stokes2cIndices) Indices;

    enum {
        lCompIdx = Indices::lCompIdx,
        gCompIdx = Indices::gCompIdx
    };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIndex) };
    enum { transportIdx = Indices::transportIdx };

public:
    /*!
     * \brief Update all additional quantities for a given control volume.
     */
    void update(const PrimaryVariables &primaryVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int vertIdx,
                bool isOldSol)
    {
        // set the mole fractions first
        completeFluidState(primaryVars, problem, element, elemGeom, vertIdx, this->fluidState(), isOldSol);

        // update vertex data for the mass and momentum balance
        ParentType::update(primaryVars,
                           problem,
                           element,
                           elemGeom,
                           vertIdx,
                           isOldSol);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(this->fluidState());

        diffCoeff_ = FluidSystem::binaryDiffusionCoefficient(this->fluidState(),
                                                             paramCache,
                                                             phaseIdx,
                                                             lCompIdx,
                                                             gCompIdx);

        Valgrind::CheckDefined(diffCoeff_);
    };

    /*!
     * \@copydoc BoxModel::completeFluidState
     */
    static void completeFluidState(const PrimaryVariables& primaryVars,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& elemGeom,
                                   int scvIdx,
                                   FluidState& fluidState,
                                   bool isOldSol = false)
    {
        Scalar massFraction[numComponents];
        massFraction[lCompIdx] = primaryVars[transportIdx];
        massFraction[gCompIdx] = 1 - massFraction[lCompIdx];

        // calculate average molar mass of the gas phase
        Scalar M1 = FluidSystem::molarMass(lCompIdx);
        Scalar M2 = FluidSystem::molarMass(gCompIdx);
        Scalar X2 = massFraction[gCompIdx];
        Scalar avgMolarMass = M1*M2/(M2 + X2*(M1 - M2));

        // convert mass to mole fractions and set the fluid state
        fluidState.setMoleFraction(phaseIdx, lCompIdx, massFraction[lCompIdx]*avgMolarMass/M1);
        fluidState.setMoleFraction(phaseIdx, gCompIdx, massFraction[gCompIdx]*avgMolarMass/M2);
    }

    /*!
     * \brief Returns the binary (mass) diffusion coefficient
     */
    Scalar diffusionCoeff() const
    { return diffCoeff_; }

protected:
    Scalar diffCoeff_; //!< Binary diffusion coefficient
};

} // end namespace

#endif
