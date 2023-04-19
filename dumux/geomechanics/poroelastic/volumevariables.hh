// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoroElastic
 * \brief Quantities required by the poroelastic model defined on a sub-control volume.
 */
#ifndef DUMUX_POROELASTIC_VOLUME_VARIABLES_HH
#define DUMUX_POROELASTIC_VOLUME_VARIABLES_HH
#include <type_traits>
#include <dune/common/exceptions.hh>


namespace Dumux {

/*!
 * \ingroup PoroElastic
 * \brief Contains the quantities which are constant within a
 *        finite volume in the poroelastic model.
 *
 * \tparam Traits Class encapsulating types to be used by the vol vars
 */
template<class Traits>
class PoroElasticVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using ModelTraits = typename Traits::ModelTraits;

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
        priVars_ = elemSol[scv.localDofIndex()];
        extrusionFactor_ = problem.spatialParams().extrusionFactor(element, scv, elemSol);

        //! set the volume fractions of the solid components
        updateSolidVolumeFractions_(elemSol, problem, element, scv);
        // set the temperature of the solid phase
        setSolidTemperature_(problem, element, scv, elemSol);
        // update the density of the solid phase
        solidState_.setDensity(SolidSystem::density(solidState_));
    }

    //! Return the average porosity \f$\mathrm{[-]}\f$ within the scv.
    Scalar solidDensity() const
    { return solidState_.density(); }

    //! Return the average porosity \f$\mathrm{[-]}\f$ within the scv
    Scalar porosity() const
    { return solidState_.porosity(); }

    //! Returns the permeability within the scv in \f$[m^2]\f$.
    Scalar displacement(unsigned int dir) const
    { return priVars_[ Indices::momentum(dir) ]; }

    //! Returns the displacement vector within the scv in \f$[m]\f$.
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
    //! updates the volume fractions of the solid components
    template<class ElemSol, class Problem, class Element, class Scv>
    void updateSolidVolumeFractions_(const ElemSol& elemSol,
                                     const Problem& problem,
                                     const Element& element,
                                     const Scv& scv)
    {
        static constexpr int numSolidComp = SolidState::numComponents;
        static constexpr int numInertComp = SolidState::numInertComponents;

        // first, set inert volume fractions from the spatial params
        const auto& sp = problem.spatialParams();
        for (int sCompIdx = numSolidComp-numInertComp; sCompIdx < numSolidComp; ++sCompIdx)
            solidState_.setVolumeFraction(sCompIdx,
                                          sp.template inertVolumeFraction<SolidSystem>(element, scv, elemSol, sCompIdx));

        // second, set the volume fractions of the (possibly) reacting components
        // these may come from a coupled flow model which considers mineralization,
        // so we make reactiveVolumeFraction a params interface in which users can
        // retrieve the current volume fractions from the flow model.
        if (!(SolidState::isInert()))
            for (int sCompIdx = 0; sCompIdx < numSolidComp-numInertComp; ++sCompIdx)
                solidState_.setVolumeFraction(sCompIdx,
                                              sp.template reactiveVolumeFraction<SolidSystem>(element, scv, elemSol, sCompIdx));
    }

    //! sets the temperature in the solid state for non-isothermal models
    static constexpr bool enableEnergyBalance = ModelTraits::enableEnergyBalance();
    template< class Problem, class Element, class Scv, class ElemSol,
              bool enableEB = enableEnergyBalance, typename std::enable_if_t<enableEB, bool> = 0 >
    void setSolidTemperature_(const Problem& problem, const Element& element, const Scv& scv, const ElemSol& elemSol)
    { DUNE_THROW(Dune::NotImplemented, "Non-isothermal poroelastic model."); }

    //! sets the temperature in the solid state for isothermal models
    template< class Problem, class Element, class Scv, class ElemSol,
              bool enableEB = enableEnergyBalance, typename std::enable_if_t<!enableEB, bool> = 0 >
    void setSolidTemperature_(const Problem& problem, const Element& element, const Scv& scv, const ElemSol& elemSol)
    { solidState_.setTemperature(problem.spatialParams().temperature(element, scv, elemSol)); }

    // data members
    Scalar extrusionFactor_;
    PrimaryVariables priVars_;
    SolidState solidState_;
};
} // end namespace Dumux

#endif
