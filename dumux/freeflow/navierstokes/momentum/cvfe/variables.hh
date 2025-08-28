// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief Variables for the single-phase Navier-Stokes model.
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_VARIABLES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_VARIABLES_HH

#include <dumux/common/concepts/ipdata_.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Variables for the single-phase Navier-Stokes model.
 */
template <class Traits>
class NavierStokesMomentumCVFEVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;

    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::dim());

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
        extrusionFactor_ = problem.spatialParams().extrusionFactor(fvGeometry, ipData, elemSol);
    }

    /*!
     * \brief Return how much the localDof is extruded.
     *
     * This means the factor by which a lower-dimensional
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    PrimaryVariables velocity() const
    { return priVars_; }

    Scalar velocity(const int dirIdx) const
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

private:
    PrimaryVariables priVars_;
    Scalar extrusionFactor_;
};

} // end namespace Dumux

#endif
