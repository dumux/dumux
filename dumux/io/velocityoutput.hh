// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
     * \brief Default constructor
     */
    VelocityOutput() = default;

    //! virtual destructor
    virtual ~VelocityOutput() {};

    //! returns whether or not velocity output is enabled
    virtual bool enableOutput() const { return false; }

    //! returns the phase name of a given phase index
    virtual std::string phaseName(int phaseIdx) const { return "none"; }

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
