// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LongWaveModel
 * \copydoc Dumux::LongWaveVolumeVariables
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_VOLUMEVARIABLES_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_VOLUMEVARIABLES_HH

namespace Dumux {

/*!
 * \ingroup LongWaveModel
 * \brief Volume variables for the long-wave approximations of the shallow water equations
 */
template <class Traits>
class LongWaveVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using Indices = typename Traits::ModelTraits::Indices;
    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());

public:
    using PrimaryVariables = typename Traits::PrimaryVariables;

    //! Update all quantities for a given control volume
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        priVars_ = elemSol[scv.indexInElement()];
        bedSurface_ = problem.spatialParams().bedSurface(element, scv);
        extrusionFactor_ = problem.spatialParams().extrusionFactor(element, scv, elemSol);
    }

    //! Return the water depth
    Scalar waterDepth() const
    { return priVars_[Indices::waterDepthIdx]; }

    //! Return the elevation of the bed surface
    Scalar bedSurface() const
    { return bedSurface_; }

    //! Return the elevation of the free surface (bed surface plus water depth)
    Scalar freeSurface() const
    { return bedSurface_ + waterDepth(); }

    //! Return a component of the primary variable vector
    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    //! Return the primary variable vector
    const PrimaryVariables& priVars() const
    { return priVars_; }

    //! Return how much the sub-control volume is extruded, e.g. the width of a one-dimensional channel
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

private:
    PrimaryVariables priVars_;
    Scalar bedSurface_;
    Scalar extrusionFactor_;
};

} // end namespace Dumux

#endif
