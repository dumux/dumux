// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 *
 * \copydoc Dumux::NavierStokesScalarConservationModelVolumeVariables
 */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2P_VOLUME_VARIABLES_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2P_VOLUME_VARIABLES_HH

#include <dumux/freeflow/navierstokes/scalarvolumevariables.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Volume variables for the multi-phase Navier-Stokes model.
 */
template <class Traits>
class NavierStokesMassTwoPVolumeVariables
: public NavierStokesScalarConservationModelVolumeVariables<Traits>
{
    using ParentType = NavierStokesScalarConservationModelVolumeVariables<Traits>;
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;
    //! export the indices type
    using Indices = typename Traits::ModelTraits::Indices;
    //! Export the underlying fluid system
    using FluidSystem = typename Traits::FluidSystem;
    using FluidState = typename Traits::FluidState;

    //! Return number of phases considered by the model
    static constexpr int numFluidPhases() { return Traits::ModelTraits::numFluidPhases(); }
    //! Return number of components considered by the model
    static constexpr int numFluidComponents() { return Traits::ModelTraits::numFluidComponents(); }

     /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        const auto& priVars = elemSol[scv.localDofIndex()];

        const auto phaseField = std::clamp<Scalar>(priVars[Indices::phaseFieldIdx], -1.0, 1.0);
        density_ = problem.mixtureDensity(phaseField);
        viscosity_ = problem.mixtureViscosity(phaseField);
    }

    /*!
     * \brief Returns the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     */
    Scalar pressure() const
    { return this->priVar(Indices::pressureIdx); }

    /*!
     * \brief Returns the phase field \f$ f \f$ within the control volume.
     */
    Scalar phaseField() const
    { return this->priVar(Indices::phaseFieldIdx); }

    /*!
     * \brief Returns the chemical potential \f$ \mu \mathrm{[J/mol]} \f$ within the control volume.
     */
    Scalar chemicalPotential() const
    { return this->priVar(Indices::chemicalPotentialIdx); }

    /*!
     * \brief Returns the dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar viscosity() const
    { return viscosity_; }

    /*!
     * \brief Returns the mass mixture density \f$\mathrm{[kg/m^3]}\f$ of a given phase within the
     *        control volume.
     */
    Scalar density() const
    { return density_; }

protected:
    Scalar density_; //!< The mass mixture density in the control volume
    Scalar viscosity_; //!< The dynamic viscosity in the control volume
};

} // end namespace Dumux

#endif
