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
 *
 * \brief Velocity output for implicit (porous media) models
 */
#ifndef DUMUX_STAGGERED_FF_VELOCITYOUTPUT_HH
#define DUMUX_STAGGERED_FF_VELOCITYOUTPUT_HH

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

namespace Properties
{
    NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output
}

/*!
 * \brief Velocity output for implicit (porous media) models
 */
template<class TypeTag>
class StaggeredFreeFlowVelocityOutput
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    static const bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using CoordScalar = typename GridView::ctype;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using ReferenceElements = Dune::ReferenceElements<CoordScalar, dim>;

public:
    /*!
     * \brief Constructor initializes the static data with the initial solution.
     *
     * \param problem The problem to be solved
     */
    StaggeredFreeFlowVelocityOutput(const Problem& problem)
    : problem_(problem)
    {
        // check, if velocity output can be used (works only for cubes so far)
        velocityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddVelocity);
    }

    bool enableOutput()
    { return velocityOutput_; }

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
        for (auto&& scv : scvs(fvGeometry))
        {
            auto dofIdxGlobal = scv.dofIndex();

            for (auto&& scvf : scvfs(fvGeometry))
            {
                auto& origFaceVars = problem_.model().curGlobalFaceVars().faceVars(scvf.dofIndex());
                auto dirIdx = scvf.directionIndex();
                velocity[dofIdxGlobal][dirIdx] += 0.5*origFaceVars.velocity();
            }
        }
    }


private:
    const Problem& problem_;
    bool velocityOutput_;
};

} // end namespace Dumux

#endif
