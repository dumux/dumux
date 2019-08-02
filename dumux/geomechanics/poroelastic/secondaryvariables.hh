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
 * \brief Quantities required by the poroelasticity
 *        model defined at an integration point.
 */
#ifndef DUMUX_POROELASTIC_SECONDARY_VARIABLES_HH
#define DUMUX_POROELASTIC_SECONDARY_VARIABLES_HH

#include <dumux/discretization/fem/secondaryvariablesbase.hh>

namespace Dumux {

/*!
 * \ingroup ElasticFemModel
 * \ingroup FemImplicitSecondaryVariables
 * \brief Contains the quantities of the linear elasticity model evaluated at an integration point.
 */
template<class Traits>
class PoroElasticSecondaryVariables : public SecondaryVariablesBase<Traits>
{
    using ParentType = SecondaryVariablesBase<Traits>;
    using Scalar = typename Traits::PrimaryVariables::value_type;

    //! The elastic model only makes sense with inert solid systems
    static_assert(Traits::SolidSystem::isInert(), "Elastic model can only be used with inert solid systems");

public:
    //! export the type used for displacement vectors
    using DisplacementVector = typename Traits::DisplacementVector;
    //! export the type encapsulating primary variable indices
    using Indices = typename Traits::ModelTraits::Indices;
    //! export type of solid state
    using SolidState = typename Traits::SolidState;
    //! export the solid system used
    using SolidSystem = typename Traits::SolidSystem;

    /*!
     * \brief Update all quantities for a given integration point
     *
     * \param elemSol The solution at the dofs connected to this element
     * \param problem The problem to be solved
     * \param element The element
     * \param ipData Container holding data on shape functions and gradients at the ip
     */
    template<class ElemSol, class Problem, class Element, class IpData>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const IpData& ipData)
    {
        ParentType::update(elemSol, problem, element, ipData);

        // update volume fractions
        updateSolidVolumeFractions_(elemSol, problem, element, ipData);
        // set the temperature of the solid phase
        setSolidTemperature_(problem);
        // update the density of the solid phase
        solidState_.setDensity(SolidSystem::density(solidState_));
    }

    //! Return the solid density \f$\mathrm{[kg/m^3]}\f$ within the control volume.
    Scalar solidDensity() const
    { return solidState_.density(); }

    //! Return the average porosity \f$\mathrm{[-]}\f$ within the scv
    Scalar porosity() const
    { return solidState_.porosity(); }

    //! Returns the displacement in a given direction in \f$[m]\f$.
    Scalar displacement(unsigned int dir) const
    { return this->priVar(Indices::momentum(dir)); }

    //! Returns the displacement vector in \f$[m]\f$.
    DisplacementVector displacement() const
    {
        DisplacementVector d;
        for (int dir = 0; dir < d.size(); ++dir)
            d[dir] = displacement(dir);
        return d;
    }

private:
    //! updates the volume fractions of the solid components
    template<class ElemSol, class Problem, class Element, class IpData>
    void updateSolidVolumeFractions_(const ElemSol& elemSol,
                                     const Problem& problem,
                                     const Element& element,
                                     const IpData& ipData)
    {
        static constexpr int numSolidComp = SolidState::numComponents;
        static constexpr int numInertComp = SolidState::numInertComponents;

        // first, set inert volume fractions from the spatial params
        const auto& sp = problem.spatialParams();
        for (int sCompIdx = numSolidComp-numInertComp; sCompIdx < numSolidComp; ++sCompIdx)
            solidState_.setVolumeFraction(sCompIdx,
                                          sp.template inertVolumeFraction<SolidSystem>(element, ipData, elemSol, sCompIdx));

        // second, set the volume fractions of the (possibly) reacting components
        if (!(SolidState::isInert()))
            for (int sCompIdx = 0; sCompIdx < numSolidComp-numInertComp; ++sCompIdx)
                solidState_.setVolumeFraction(sCompIdx,
                                              sp.template reactiveVolumeFraction<SolidSystem>(element, ipData, elemSol, sCompIdx));
    }

    //! sets the temperature in the solid state for non-isothermal models
    static constexpr bool enableEnergyBalance = Traits::ModelTraits::enableEnergyBalance();
    template< class Problem,
              bool enableEB = enableEnergyBalance, typename std::enable_if_t<enableEB, bool> = 0 >
    void setSolidTemperature_(const Problem& problem)
    { DUNE_THROW(Dune::InvalidStateException, "Non-isothermal elastic model."); }

    //! sets the temperature in the solid state for isothermal models
    template< class Problem,
              bool enableEB = enableEnergyBalance, typename std::enable_if_t<!enableEB, bool> = 0 >
    void setSolidTemperature_(const Problem& problem)
    { solidState_.setTemperature(problem.temperature()); }

    SolidState solidState_;
};

}

#endif
