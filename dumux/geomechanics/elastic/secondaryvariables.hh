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
 * \brief Quantities required by the linear elasticity box
 *        model defined on a vertex.
 */
#ifndef DUMUX_ELASTIC_SECONDARY_VARIABLES_HH
#define DUMUX_ELASTIC_SECONDARY_VARIABLES_HH

namespace Dumux {

/*!
 * \ingroup ElasticFemModel
 * \ingroup FemImplicitSecondaryVariables
 * \brief Contains the quantities of the linear elasticity model evaluated at an integration point.
 */
template<class Traits>
class ElasticSecondaryVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using ModelTraits = typename Traits::ModelTraits;

    //! The elastic model only makes sense with inert solid systems
    static_assert(Traits::SolidSystem::isInert(), "Elastic model can only be used with inert solid systems");

public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;
    //! export the type used for displacement vectors
    using DisplacementVector = typename Traits::DisplacementVector;
    //! export the type encapsulating primary variable indices
    using Indices = typename ModelTraits::Indices;
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
        // interpolate primary variables
        priVars_ = 0.0;
        for (unsigned int i = 0; i < elemSol.size(); ++i)
        {
            PrimaryVariables tmp(elemSol[i]);
            tmp *= ipData.shapeValues()[i];
            priVars_ += tmp;
        }

        extrusionFactor_ = problem.extrusionFactor(element, ipData, elemSol);

        //! set the volume fractions of the solid components
        for (int sCompIdx = solidState_.numComponents-solidState_.numInertComponents; sCompIdx < solidState_.numComponents; ++sCompIdx)
            solidState_.setVolumeFraction(sCompIdx,
                                          problem.spatialParams().template inertVolumeFraction<SolidSystem>(element, ipData, elemSol, sCompIdx));

        // TODO: what to do for non-inert solids?
        // if (!(solidState.isInert())) {}

        // set the temperature of the solid phase
        setSolidTemperature_(problem);
        // update the density of the solid phase
        solidState_.setDensity(SolidSystem::density(solidState_));
    }

    //! Return the solid density \f$\mathrm{[kg/m^3]}\f$ within the control volume.
    Scalar solidDensity() const
    { return solidState_.density(); }

    //! Returns the displacement in a given direction in \f$[m]\f$.
    Scalar displacement(unsigned int dir) const
    { return priVars_[ Indices::momentum(dir) ]; }

    //! Returns the displacement vector in \f$[m]\f$.
    DisplacementVector displacement() const
    {
        DisplacementVector d;
        for (int dir = 0; dir < d.size(); ++dir)
            d[dir] = displacement(dir);
        return d;
    }

    //! Return a component of primary variable vector for a given index
    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    //! Return the vector of primary variables
    const PrimaryVariables& priVars() const
    { return priVars_; }

    //! TODO We don't know yet how to interpret extrusion for mechanics
    static constexpr Scalar extrusionFactor()
    { return 1.0; }

private:
    //! sets the temperature in the solid state for non-isothermal models
    static constexpr bool enableEnergyBalance = ModelTraits::enableEnergyBalance();
    template< class Problem,
              bool enableEB = enableEnergyBalance, typename std::enable_if_t<enableEB, bool> = 0 >
    void setSolidTemperature_(const Problem& problem)
    { DUNE_THROW(Dune::InvalidStateException, "Non-isothermal elastic model."); }

    //! sets the temperature in the solid state for isothermal models
    template< class Problem,
              bool enableEB = enableEnergyBalance, typename std::enable_if_t<!enableEB, bool> = 0 >
    void setSolidTemperature_(const Problem& problem)
    { solidState_.setTemperature(problem.temperature()); }

    // data members
    Scalar extrusionFactor_;
    PrimaryVariables priVars_;
    SolidState solidState_;
};

}

#endif
