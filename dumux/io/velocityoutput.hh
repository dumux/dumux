// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Default velocity output policy for porous media models
 */
#ifndef DUMUX_IO_VELOCITYOUTPUT_HH
#define DUMUX_IO_VELOCITYOUTPUT_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Velocity output for implicit (porous media) models
 */
template<class GridVariables>
class VelocityOutput
{
    using Scalar = typename GridVariables::Scalar;
    static constexpr int dimWorld = GridVariables::GridGeometry::GridView::dimensionworld;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVarsCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using FVElementGeometry = typename GridVariables::GridGeometry::LocalView;
    using Element = typename GridVariables::GridGeometry::GridView::template Codim<0>::Entity;

public:
    using VelocityVector = std::vector<Dune::FieldVector<Scalar, dimWorld>>;

    /*!
    * \brief A container for possible velocity data types
    */
    enum class FieldType
    {
        element, vertex, automatic
    };

    /*!
     * \brief Default constructor
     */
    VelocityOutput() = default;

    //! virtual destructor
    virtual ~VelocityOutput() {};

    //! returns whether or not velocity output is enabled
    virtual bool enableOutput() const { return false; }

    //! returns the phase name of a given phase index
    virtual std::string phaseName(int phaseIdx) const { return "none"; }

    //! returns the field type
    virtual FieldType fieldType() const { return FieldType::automatic; }

    //! returns the number of phases
    virtual int numFluidPhases() const { return 0; }

    //! Calculate the velocities for the scvs in the element
    //! We assume the local containers to be bound to the complete stencil
    virtual void calculateVelocity(VelocityVector& velocity,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars,
                                   const ElementFluxVarsCache& elemFluxVarsCache,
                                   int phaseIdx) const
    {}
};

} // end namespace Dumux

#endif
