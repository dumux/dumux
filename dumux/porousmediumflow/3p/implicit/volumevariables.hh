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
 * \brief Contains the quantities which are constant within a
 *        finite volume in the three-phase model.
 */
#ifndef DUMUX_3P_VOLUME_VARIABLES_HH
#define DUMUX_3P_VOLUME_VARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

/*!
 * \ingroup ThreePModel
 * \ingroup PorousMediumFlowVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in three-phase model.
 */
template <class TypeTag>
class ThreePVolumeVariables : public PorousMediumFlowVolumeVariables<TypeTag>
{
    using ParentType = PorousMediumFlowVolumeVariables<TypeTag>;

    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using PermeabilityType = typename SpatialParams::PermeabilityType;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    enum {
        dim = GridView::dimension,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),

        wPhaseIdx = Indices::wPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        swIdx = Indices::swIdx,
        snIdx = Indices::snIdx,
        pressureIdx = Indices::pressureIdx
    };

    static const Scalar R; // universal gas constant

    enum { isBox = GET_PROP_VALUE(TypeTag, DiscretizationMethod) == DiscretizationMethods::Box };
    enum { dofCodim = isBox ? dim : 0 };

public:
    //! The type of the object returned by the fluidState() method
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;


    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        // capillary pressure parameters
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        completeFluidState(elemSol, problem, element, scv, fluidState_);

         // mobilities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            mobility_[phaseIdx] = MaterialLaw::kr(materialParams, phaseIdx,
                                 fluidState_.saturation(wPhaseIdx),
                                 fluidState_.saturation(nPhaseIdx),
                                 fluidState_.saturation(gPhaseIdx))
                                 / fluidState_.viscosity(phaseIdx);
            Valgrind::CheckDefined(mobility_[phaseIdx]);
        }

        // porosity
        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);

        Valgrind::CheckDefined(porosity_);
        Valgrind::CheckDefined(permeability_);
    }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     */
    static void completeFluidState(const ElementSolutionVector& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        Scalar t = ParentType::temperature(elemSol, problem, element, scv);
        fluidState.setTemperature(t);

        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv);

        const Scalar sw = priVars[swIdx];
        const Scalar sn = priVars[snIdx];
        const Scalar sg = 1.0 - sw - sn;

        Valgrind::CheckDefined(sg);

        fluidState.setSaturation(wPhaseIdx, sw);
        fluidState.setSaturation(gPhaseIdx, sg);
        fluidState.setSaturation(nPhaseIdx, sn);

        /* now the pressures */
        const Scalar pg = priVars[pressureIdx];

        // calculate capillary pressures
        const Scalar pcgw = MaterialLaw::pcgw(materialParams, sw);
        const Scalar pcnw = MaterialLaw::pcnw(materialParams, sw);
        const Scalar pcgn = MaterialLaw::pcgn(materialParams, sw + sn);

        const Scalar pcAlpha = MaterialLaw::pcAlpha(materialParams, sn);
        const Scalar pcNW1 = 0.0; // TODO: this should be possible to assign in the problem file

        const Scalar pn = pg- pcAlpha * pcgn - (1.0 - pcAlpha)*(pcgw - pcNW1);
        const Scalar pw = pn - pcAlpha * pcnw - (1.0 - pcAlpha)*pcNW1;

        fluidState.setPressure(wPhaseIdx, pw);
        fluidState.setPressure(gPhaseIdx, pg);
        fluidState.setPressure(nPhaseIdx, pn);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // compute and set the viscosity
            const Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx,mu);

            // compute and set the density
            const Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);

            // compute and set the enthalpy
            const Scalar h = ParentType::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the effective saturation of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(const int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperatures of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(const int phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \brief Returns the effective capillary pressure within the control volume.
     */
    Scalar capillaryPressure() const
    { return fluidState_.capillaryPressure(); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

protected:

    Scalar porosity_;        //!< Effective porosity within the control volume
    PermeabilityType permeability_;
    Scalar mobility_[numPhases];  //!< Effective mobility within the control volume
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

template <class TypeTag>
const typename ThreePVolumeVariables<TypeTag>::Scalar ThreePVolumeVariables<TypeTag>::R = Constants<typename GET_PROP_TYPE(TypeTag, Scalar)>::R;

} // end namespace

#endif
