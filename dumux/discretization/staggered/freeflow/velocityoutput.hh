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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFreeFlowVelocityOutput
 */
#ifndef DUMUX_STAGGERED_FF_VELOCITYOUTPUT_HH
#define DUMUX_STAGGERED_FF_VELOCITYOUTPUT_HH

#include <dune/common/fvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

namespace Dumux
{

/*!
 * \ingroup StaggeredDiscretization
 * \brief Velocity output for staggered free-flow models
 */
template<class TypeTag>
class StaggeredFreeFlowVelocityOutput
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

    using Element = typename GridView::template Codim<0>::Entity;
    using CoordScalar = typename GridView::ctype;

public:
    /*!
     * \brief Constructor initializes the static data with the initial solution.
     *
     * \param problem The problem
     * \param fvGridGeometry The fvGridGeometry
     * \param gridVariables The gridVariables
     * \param sol The solution vector
     */
    StaggeredFreeFlowVelocityOutput(const Problem& problem,
                           const FVGridGeometry& fvGridGeometry,
                           const GridVariables& gridVariables,
                           const SolutionVector& sol)
    : problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    , gridVariables_(gridVariables)
    , sol_(sol)
    {
        // check if velocity vectors shall be written to the VTK file
        velocityOutput_ = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Vtk.AddVelocity");
    }

    //! Returns whether to enable the velocity output or not
    bool enableOutput()
    { return velocityOutput_; }

    // returns the name of the phase for a given index
    static std::string phaseName(int phaseIdx)
    { return GET_PROP_TYPE(TypeTag, FluidSystem)::phaseName(phaseIdx); }

    //! Return the problem boundary types
    auto problemBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    { return problem_.boundaryTypes(element, scvf); }

    //! Calculate the velocities for the scvs in the element
    //! We assume the local containers to be bound to the complete stencil
    template<class VelocityVector>
    void calculateVelocity(VelocityVector& velocity,
                           const ElementVolumeVariables& elemVolVars,
                           const FVElementGeometry& fvGeometry,
                           const Element& element,
                           int phaseIdx)
    {
        auto elemFaceVars = localView(gridVariables_.curGridFaceVars());
        elemFaceVars.bindElement(element, fvGeometry, sol_);
        for (auto&& scv : scvs(fvGeometry))
        {
            auto dofIdxGlobal = scv.dofIndex();

            for (auto&& scvf : scvfs(fvGeometry))
            {
                auto dirIdx = scvf.directionIndex();
                velocity[dofIdxGlobal][dirIdx] += 0.5*elemFaceVars[scvf].velocitySelf();
            }
        }
    }

private:
    const Problem& problem_;
    const FVGridGeometry& fvGridGeometry_;
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;
    bool velocityOutput_;
};

} // end namespace Dumux

#endif
