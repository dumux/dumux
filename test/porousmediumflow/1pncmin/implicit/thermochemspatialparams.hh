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
 * \ingroup OnePNCMinTests
 * \brief Definition of the spatial parameters for the thermochemistry
 *        problem which uses the non-insothermal 1pncmin model
 */

#ifndef DUMUX_THERMOCHEM_SPATIAL_PARAMS_HH
#define DUMUX_THERMOCHEM_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

#include <dumux/material/fluidmatrixinteractions/porosityreactivebed.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>
#include <dumux/material/fluidmatrixinteractions/mineralization/effectivesoliddensity.hh>
#include <dumux/material/fluidmatrixinteractions/mineralization/effectivesolidheatcapacity.hh>

namespace Dumux {

//forward declaration
template<class TypeTag>
class ThermoChemSpatialParams;

namespace Properties {

// The spatial parameters TypeTag
NEW_TYPE_TAG(ThermoChemSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(ThermoChemSpatialParams, SpatialParams, Dumux::ThermoChemSpatialParams<TypeTag>);

} // end namespace Properties

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup BoxTestProblems
 * \brief Definition of the spatial parameters for the FuelCell
 *        problem which uses the isothermal 2p2c box model
 */
template<class TypeTag>
class ThermoChemSpatialParams : public FVSpatialParamsOneP<TypeTag>
{
    using ParentType = FVSpatialParamsOneP<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    enum {
        dimWorld=GridView::dimensionworld,

        numSPhases =  ModelTraits::numSPhases(),
        numPhases = ModelTraits::numPhases(),
        cPhaseIdx = FluidSystem::cPhaseIdx,
        hPhaseIdx = FluidSystem::hPhaseIdx
    };

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    using PorosityLaw = PorosityReactiveBed<TypeTag>;
    using PermeabilityLaw = PermeabilityKozenyCarman<TypeTag>;
    using EffectiveSolidRho = EffectiveSolidDensity<TypeTag>;
    using EffectiveSolidCp = EffectiveSolidHeatCapacity<TypeTag>;

public:
    // type used for the permeability (i.e. tensor or scalar)
    using PermeabilityType = Scalar;
    /*!
     * \brief The constructor
     *
     * \param problem The Problem
     */
    ThermoChemSpatialParams(const Problem& problem)
    : ParentType(problem)
    {
        //thermal conductivity of CaO
        lambdaSolid_ = 0.4; //[W/(m*K)] Nagel et al [2013b]
        rho_[cPhaseIdx-numPhases] = 1656; //[kg/m^3] density of CaO (see Shao et al. 2014)
        rho_[hPhaseIdx-numPhases] = 2200; //[kg/m^3] density of Ca(OH)_2 (see Shao et al. 2014)
        cp_[cPhaseIdx-numPhases] = 934; //[J/kgK] heat capacity of CaO (see Shao et al. 2014)
        cp_[hPhaseIdx-numPhases] = 1530; //[J/kgK] heat capacity of Ca(OH)_2 (see Shao et al. 2014)
        eps_ = 1e-6;
        poroLaw_.init(*this);
        permLaw_.init(*this);
        effSolRho_.init(*this);
        effSolCp_.init(*this);
    }

    /*!
     *  \brief Define the reference permeability \f$[m^2]\f$ distribution
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar referencePermeability(const Element& element, const SubControlVolume &scv) const
    { return 8.53e-12; }

    /*!
     *  \brief Define the reference porosity \f$[-]\f$ distribution
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar referencePorosity(const Element& element, const SubControlVolume &scv) const
    {
        return  0.8;
    }

    /*! Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param element The finite volume element
     *  \param scv The sub-control volume
     *  \param elemSol The element solution
     *
     *  Solution dependent permeability function
     */
    template<class ElementSolution>
    Scalar permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    { return permLaw_.evaluatePermeability(element, scv, elemSol); }

    /*!
     *  \brief Define the minimum porosity \f$[-]\f$ distribution
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     */
    Scalar minPorosity(const Element& element, const SubControlVolume &scv) const
    {
//         return 0.604; //intrinsic porosity of CaO2H2 (see Nagel et al. 2014)
        return 0.8;
    }

    /*!
     *  \brief Define the minimum porosity \f$[-]\f$ after clogging caused by mineralization
     *
     *  \param element The finite element
     *  \param scv The sub-control volume
     *  \param elemSol The element solution
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
//         return poroLaw_.evaluatePorosity(element, scv, elemSol);
        return 0.8;
    }

    /*!
     * \brief Returns the average heat capacity \f$[J / (kg K)]\f$ of solid phases.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    template<class ElementSolution>
    Scalar solidHeatCapacity(const Element &element,
                             const SubControlVolume& scv,
                             const ElementSolution& elemSol) const
    {
        return effSolCp_.effectiveSolidHeatCapacity(element, scv, elemSol);
    }

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the pure solid phases.
     */
    template<class ElementSolution>
    Scalar solidPhaseHeatCapacity(const Element &element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol,
                                  int sPhaseIdx) const
    {
        return cp_[sPhaseIdx];
    }

    /*!
     * \brief Returns the average mass density \f$[kg / m^3]\f$ of the solid phases.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    template<class ElementSolution>
    Scalar solidDensity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    {
        return effSolRho_.effectiveSolidDensity(element, scv, elemSol);
    }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the pure solid phases.
     */
    template<class ElementSolution>
    Scalar solidPhaseDensity(const Element &element,
                             const SubControlVolume& scv,
                             const ElementSolution& elemSol,
                             int sPhaseIdx) const
    {
        return rho_[sPhaseIdx];
    }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     */
    template<class ElementSolution>
    Scalar solidThermalConductivity(const Element &element,
                                    const SubControlVolume& scv,
                                    const ElementSolution& elemSol) const
    { return lambdaSolid_; }

private:

   Scalar eps_;
   Scalar lambdaSolid_;
   std::array<Scalar, numSPhases> rho_;
   std::array<Scalar, numSPhases> cp_;
   PorosityLaw poroLaw_;
   PermeabilityLaw permLaw_;
   EffectiveSolidRho effSolRho_;
   EffectiveSolidCp effSolCp_;
};

}//end namespace

#endif
