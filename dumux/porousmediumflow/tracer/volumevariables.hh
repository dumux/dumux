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
 * \ingroup TracerModel
 * \brief Quantities required by the tracer model in a control volume
 */
#ifndef DUMUX_TRACER_VOLUME_VARIABLES_HH
#define DUMUX_TRACER_VOLUME_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/volumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup TracerModel
 * \brief Contains the quantities which are constant within a
 *        finite volume for the tracer model.
 */
template <class TypeTag>
class TracerVolumeVariables : public PorousMediumFlowVolumeVariables<TypeTag>
{
    using ParentType = PorousMediumFlowVolumeVariables<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

public:

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    {
        // update parent type sets primary variables
        ParentType::update(elemSol, problem, element, scv);

        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
        dispersivity_ = problem.spatialParams().dispersivity(element, scv, elemSol);

        // the spatial params special to the tracer model
        fluidDensity_ = problem.spatialParams().fluidDensity(element, scv);
        fluidMolarMass_ = problem.spatialParams().fluidMolarMass(element, scv);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            moleOrMassFraction_[compIdx] = this->priVars()[compIdx];
            diffCoeff_[compIdx] = FluidSystem::binaryDiffusionCoefficient(compIdx, problem, element, scv);
        }
    }

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     *
     * \param pIdx TODO docme!
     */
    Scalar density(int pIdx = 0) const
    { return fluidDensity_; }

    /*!
     * \brief Return the saturation
     *
     * This method is here for compatibility reasons with other models. The saturation
     * is always 1.0 in a one-phasic context.
     *
     * \param pIdx TODO docme!
     */
    Scalar saturation(int pIdx = 0) const
    { return 1.0; }

    /*!
     * \brief Return the mobility
     *
     * This method is here for compatibility reasons with other models. The mobility is always 1
     * for one-phasic models where the velocity field is given
     *
     * \param pIdx TODO docme!
     */
    Scalar mobility(int pIdx = 0) const
    { return 1.0; }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * \param pIdx TODO docme!
     */
    Scalar molarDensity(int pIdx = 0) const
    { return fluidDensity_/fluidMolarMass_; }

    /*!
     * \brief Return mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
     *
     * \param pIdx TODO docme!
     * \param compIdx The index of the component
     */
    Scalar moleFraction(int pIdx, int compIdx) const
    { return useMoles ? moleOrMassFraction_[compIdx] : moleOrMassFraction_[compIdx]/FluidSystem::molarMass(compIdx)*fluidMolarMass_; }

    /*!
     * \brief Return mass fraction \f$\mathrm{[kg/kg]}\f$ of a component in the phase.
     *
     * \param pIdx TODO docme!
     * \param compIdx The index of the component
     */
    Scalar massFraction(int pIdx, int compIdx) const
    { return useMoles ? moleOrMassFraction_[compIdx]*FluidSystem::molarMass(compIdx)/fluidMolarMass_ : moleOrMassFraction_[compIdx]; }

    /*!
     * \brief Return concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param pIdx TODO docme!
     * \param compIdx The index of the component
     */
    Scalar molarity(int pIdx, int compIdx) const
    { return moleFraction(pIdx, compIdx)*molarDensity(); }

    /*!
     * \brief Return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ in the fluid.
     *
     * \param pIdx TODO docme!
     * \param compIdx The index of the component
     */
    Scalar diffusionCoefficient(int pIdx, int compIdx) const
    { return diffCoeff_[compIdx]; }

    /*!
     * \brief Returns the dispersivity of the fluid's streamlines.
     */
    const GlobalPosition &dispersivity() const
    { return dispersivity_; }

    /*!
     * \brief Return the average porosity \f$\mathrm{[-]}\f$ within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

protected:
    Scalar porosity_;    // Effective porosity within the control volume
    Scalar fluidDensity_, fluidMolarMass_;
    GlobalPosition dispersivity_;
    std::array<Scalar, numComponents> diffCoeff_;
    std::array<Scalar, numComponents> moleOrMassFraction_;
};

} // end namespace Dumux

#endif
