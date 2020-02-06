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
 * \ingroup SolidEnergyModel
 * \brief Class for computation of all volume averaged quantities
 */
#ifndef DUMUX_SOLID_ENERGY_VOLUME_VARIABLES_HH
#define DUMUX_SOLID_ENERGY_VOLUME_VARIABLES_HH

#include <type_traits>

#include <dumux/material/solidstates/updatesolidvolumefractions.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup SolidEnergyModel
 * \brief Class for computation of all volume averaged quantities
 */
template<class Traits>
class SolidEnergyVolumeVariables : public PorousMediumFlowVolumeVariables<Traits>
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    static constexpr int temperatureIdx = Traits::ModelTraits::Indices::temperatureIdx;

public:
    //! export the type used for the solid state
    using SolidState = typename Traits::SolidState;
    //! export the type used for the solid system
    using SolidSystem = typename Traits::SolidSystem;

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);
        updateTemperature(elemSol, problem, element, scv, solidState_);
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, 0);
        updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
    }

    //! Fill temperature in the solid state
    template<class ElemSol, class Problem, class Element, class Scv>
    void updateTemperature(const ElemSol& elemSol,
                           const Problem& problem,
                           const Element& element,
                           const Scv& scv,
                           SolidState& solidState)
    {
        const Scalar T = elemSol[scv.localDofIndex()][temperatureIdx];
        solidState.setTemperature(T);
    }

    //! Fill solid matrix parameters in the solid state
    template<class ElemSol, class Problem, class Element, class Scv>
    void updateSolidEnergyParams(const ElemSol &elemSol,
                                 const Problem& problem,
                                 const Element &element,
                                 const Scv &scv,
                                 SolidState & solidState)
    {
        Scalar cs = solidHeatCapacity_(elemSol, problem, element, scv, solidState);
        solidState.setHeatCapacity(cs);

        Scalar rhos = solidDensity_(elemSol, problem, element, scv, solidState);
        solidState.setDensity(rhos);

        Scalar lambdas = solidThermalConductivity_(elemSol, problem, element, scv, solidState);
        solidState.setThermalConductivity(lambdas);
    }

    /*!
     * \brief Returns the temperature in the sub-control volume.
     */
    Scalar temperatureSolid() const
    { return solidState_.temperature(); }

    /*!
     * \brief Returns the temperature in the sub-control volume.
     */
    Scalar temperature() const
    { return solidState_.temperature(); }

    /*!
     * \brief Returns the total heat capacity \f$\mathrm{[J/(kg K)]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar solidHeatCapacity() const
    { return solidState_.heatCapacity(); }

    /*!
     * \brief Returns the mass density \f$\mathrm{[kg/m^3]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar solidDensity() const
    {  return  solidState_.density(); }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ of the solid phase in the sub-control volume.
     */
    Scalar solidThermalConductivity() const
    { return solidState_.thermalConductivity(); }

    /*!
     * \brief Returns the effective thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
     *        of the solid phase in the sub-control volume. In this case (non-porous) identical to the solidThermalCondutivity.
     */
    Scalar effectiveThermalConductivity() const
    { return solidThermalConductivity(); }

    /*!
     * \brief Return the average porosity \f$\mathrm{[-]}\f$ within the control volume.
     */
    Scalar porosity() const
    { return solidState_.porosity(); }

private:
    //! the solid state
    SolidState solidState_;

    /*!
     * It has to be decided if the full solid system / solid state interface is used (general option, but more complicated),
     * or the simple nonisothermal spatial params interface (simpler but less general).
     * In the simple nonisothermal spatial params interface the functions solidHeatCapacity, solidDensity, and solidThermalConductivity
     * in the spatial params overwrite the parameters given in the solid system. This only makes sense in combination
     * with the simplest solid system InertSolidPhase, and can be used to quickly change parameters in certain domain regions.
     * For setups with more general solids with several components these functions should not exist. Instead, the solid system
     * determines the values for solidHeatCapacity, solidDensity, and solidThermalConductivity depending on the given composition.
     */

    /*!
     * \name Access functions for the solidsystem / solidstate interface
     */
    // \{

    /*!
     * \brief get the solid heat capacity in an scv
     * \param elemSol the element solution vector
     * \param problem the problem to solve
     * \param element the element (codim-0-entity) the scv belongs to
     * \param scv the sub control volume
     * \param solidState the solid state
     * \note this gets selected if the user uses the solidsystem / solidstate interface
     */
    template<class ElemSol, class Problem, class Element, class Scv,
             std::enable_if_t<!Detail::hasSolidHeatCapacity<typename Problem::SpatialParams, Element, Scv, ElemSol, SolidState>(), int> = 0>
    Scalar solidHeatCapacity_(const ElemSol& elemSol,
                              const Problem& problem,
                              const Element& element,
                              const Scv& scv,
                              const SolidState& solidState)
    {
        return SolidSystem::heatCapacity(solidState);
    }

    /*!
     * \brief get the solid density in an scv
     * \param elemSol the element solution vector
     * \param problem the problem to solve
     * \param element the element (codim-0-entity) the scv belongs to
     * \param scv the sub control volume
     * \param solidState the solid state
     * \note this gets selected if the user uses the solidsystem / solidstate interface
     */
    template<class ElemSol, class Problem, class Element, class Scv,
             std::enable_if_t<!Detail::hasSolidDensity<typename Problem::SpatialParams, Element, Scv, ElemSol, SolidState>(), int> = 0>
    Scalar solidDensity_(const ElemSol& elemSol,
                         const Problem& problem,
                         const Element& element,
                         const Scv& scv,
                         const SolidState& solidState)
    {
        return SolidSystem::density(solidState);
    }

    /*!
     * \brief get the solid's thermal conductivity in an scv
     * \param elemSol the element solution vector
     * \param problem the problem to solve
     * \param element the element (codim-0-entity) the scv belongs to
     * \param scv the sub control volume
     * \param solidState the solid state
     * \note this gets selected if the user uses the solidsystem / solidstate interface
     */
    template<class ElemSol, class Problem, class Element, class Scv,
             std::enable_if_t<!Detail::hasSolidThermalConductivity<typename Problem::SpatialParams, Element, Scv, ElemSol, SolidState>(), int> = 0>
    Scalar solidThermalConductivity_(const ElemSol& elemSol,
                                     const Problem& problem,
                                     const Element& element,
                                     const Scv& scv,
                                     const SolidState& solidState)
    {
        return SolidSystem::thermalConductivity(solidState);
    }

    // \}

    /*!
     * \name Access functions for the simple nonisothermal spatial params interface in
     *       combination with an InertSolidPhase as solid system
     */
    // \{

    /*!
     * \brief get the solid heat capacity in an scv
     * \param elemSol the element solution vector
     * \param problem the problem to solve
     * \param element the element (codim-0-entity) the scv belongs to
     * \param scv the sub control volume
     * \param solidState the solid state
     * \note this gets selected if the user uses the simple spatial params interface in
     *       combination with an InertSolidPhase as solid system
     */
    template<class ElemSol, class Problem, class Element, class Scv,
             std::enable_if_t<Detail::hasSolidHeatCapacity<typename Problem::SpatialParams, Element, Scv, ElemSol, SolidState>(), int> = 0>
    Scalar solidHeatCapacity_(const ElemSol& elemSol,
                              const Problem& problem,
                              const Element& element,
                              const Scv& scv,
                              const SolidState& solidState)
    {
        static_assert(Detail::isInertSolidPhase<SolidSystem>::value,
            "solidHeatCapacity can only be overwritten in the spatial params when the solid system is a simple InertSolidPhase\n"
            "If you select a proper solid system, the solid heat capacity will be computed as stated in the solid system!");
        return problem.spatialParams().solidHeatCapacity(element, scv, elemSol, solidState);
    }

    /*!
     * \brief get the solid density in an scv
     * \param elemSol the element solution vector
     * \param problem the problem to solve
     * \param element the element (codim-0-entity) the scv belongs to
     * \param scv the sub control volume
     * \param solidState the solid state
     * \note this gets selected if the user uses the simple spatial params interface in
     *       combination with an InertSolidPhase as solid system
     */
    template<class ElemSol, class Problem, class Element, class Scv,
             std::enable_if_t<Detail::hasSolidDensity<typename Problem::SpatialParams, Element, Scv, ElemSol, SolidState>(), int> = 0>
    Scalar solidDensity_(const ElemSol& elemSol,
                         const Problem& problem,
                         const Element& element,
                         const Scv& scv,
                         const SolidState& solidState)
    {
        static_assert(Detail::isInertSolidPhase<SolidSystem>::value,
            "solidDensity can only be overwritten in the spatial params when the solid system is a simple InertSolidPhase\n"
            "If you select a proper solid system, the solid density will be computed as stated in the solid system!");
        return problem.spatialParams().solidDensity(element, scv, elemSol, solidState);
    }

    /*!
     * \brief get the solid's heat capacity in an scv
     * \param elemSol the element solution vector
     * \param problem the problem to solve
     * \param element the element (codim-0-entity) the scv belongs to
     * \param scv the sub control volume
     * \param solidState the solid state
     * \note this gets selected if the user uses the simple spatial params interface in
     *       combination with an InertSolidPhase as solid system
     */
    template<class ElemSol, class Problem, class Element, class Scv,
             std::enable_if_t<Detail::hasSolidThermalConductivity<typename Problem::SpatialParams, Element, Scv, ElemSol, SolidState>(), int> = 0>
    Scalar solidThermalConductivity_(const ElemSol& elemSol,
                                     const Problem& problem,
                                     const Element& element,
                                     const Scv& scv,
                                     const SolidState& solidState)
    {
        static_assert(Detail::isInertSolidPhase<SolidSystem>::value,
            "solidThermalConductivity can only be overwritten in the spatial params when the solid system is a simple InertSolidPhase\n"
            "If you select a proper solid system, the solid thermal conductivity will be computed as stated in the solid system!");
        return problem.spatialParams().solidThermalConductivity(element, scv, elemSol, solidState);
    }

    // \}

};

} // end namespace Dumux

#endif
