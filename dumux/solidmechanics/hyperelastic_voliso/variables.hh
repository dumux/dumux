// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HyperelasticVolIso
 * \brief Volume variables for the volumetric-isochoric split hyperelastic mixed u–p model.
 */
#ifndef DUMUX_SOLIDMECHANICS_HYPERELASTIC_VOLISO_VARIABLES_HH
#define DUMUX_SOLIDMECHANICS_HYPERELASTIC_VOLISO_VARIABLES_HH

#include <dumux/common/concepts/ipdata_.hh>

namespace Dumux {

/*!
 * \ingroup HyperelasticVolIso
 * \brief Volume variables for the momentum subdomain (displacement).
 */
template <class Traits>
class HyperelasticVolIsoMomentumVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;

    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());

public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;

    //! export the indices type
    using Indices = typename Traits::ModelTraits::Indices;

    /*!
     * \brief Update all quantities for a local dof
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param fvGeometry The local geometry
     * \param ipData The interpolation point data
     */
    template<class ElementSolution, class Problem, class FVElementGeometry, Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const FVElementGeometry& fvGeometry,
                const IpData& ipData)
    {
        priVars_ = elemSol[ipData.localDofIndex()];
    }

    PrimaryVariables displacement() const
    { return priVars_; }

    Scalar displacement(const int dirIdx) const
    { return priVars_[dirIdx]; }

    /*!
     * \brief Return a component of primary variable vector
     * \param pvIdx The index of the primary variable of interest
     */
    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    /*!
     * \brief Return the primary variable vector
     */
    const PrimaryVariables& priVars() const
    { return priVars_; }

    /*!
     * \brief Returns the factor by which a lower-dimensional entity needs to be extruded.
     */
    Scalar extrusionFactor() const
    { return 1.0; }

private:
    PrimaryVariables priVars_;
};


/*!
 * \ingroup HyperelasticVolIso
 * \brief Volume variables for the pressure subdomain (bulk pressure).
 */
template <class Traits>
class HyperelasticVolIsoPressureVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;

    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());

public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;

    //! export the indices type
    using Indices = typename Traits::ModelTraits::Indices;

    /*!
     * \brief Update all quantities for a local dof
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param fvGeometry The local geometry
     * \param ipData The interpolation point data
     */
    template<class ElementSolution, class Problem, class FVElementGeometry, Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const FVElementGeometry& fvGeometry,
                const IpData& ipData)
    {
        priVars_ = elemSol[ipData.localDofIndex()];
    }

    Scalar pressure() const
    { return priVars_[0]; }

    /*!
     * \brief Return a component of primary variable vector
     * \param pvIdx The index of the primary variable of interest
     */
    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    /*!
     * \brief Return the primary variable vector
     */
    const PrimaryVariables& priVars() const
    { return priVars_; }

    /*!
     * \brief Returns the factor by which a lower-dimensional entity needs to be extruded.
     */
    Scalar extrusionFactor() const
    { return 1.0; }

private:
    PrimaryVariables priVars_;
};

} // end namespace Dumux

#endif
