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
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
#ifndef DUMUX_1P_TEST_SPATIALPARAMS_HH
#define DUMUX_1P_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>

namespace Dumux
{

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
template<class TypeTag>
class OnePTestSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    OnePTestSpatialParams(const Problem& problem, const GridView& gridView)
        : ParentType(problem, gridView)
    {
        eps_ = 1.5e-7;

            lensLowerLeft_[0] = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.LensLowerLeftX);
            if (dimWorld > 1)
                lensLowerLeft_[1] = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.LensLowerLeftY);
            if (dimWorld > 2)
                lensLowerLeft_[2] = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.LensLowerLeftZ);

            lensUpperRight_[0] = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.LensUpperRightX);
            if (dimWorld > 1)
                lensUpperRight_[1] = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.LensUpperRightY);
            if (dimWorld > 2)
                lensUpperRight_[2] = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.LensUpperRightZ);

            permeability_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Permeability);
            permeabilityLens_=GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PermeabilityLens);
    }

    /*!
     * \brief Return the intrinsic permeability for the current sub-control volume in [m^2].
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index sub-control volume face where the
     *                      intrinsic velocity ought to be calculated.
     */
    Scalar intrinsicPermeability(const SubControlVolume& scv,
                                 const VolumeVariables& volVars = VolumeVariables()) const
    {
        if (isInLens_(scv.dofPosition()))
            return permeabilityLens_;
        else
            return permeability_;
    }

    /*! \brief Define the porosity in [-].
   *
   * \param element The finite element
   * \param fvGeometry The finite volume geometry
   * \param scvIdx The local index of the sub-control volume where
   */
    Scalar porosity(const SubControlVolume &scv) const
    { return 0.4; }

private:
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i) {
            if (globalPos[i] < lensLowerLeft_[i] + eps_ || globalPos[i] > lensUpperRight_[i] - eps_)
                return false;
        }
        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar permeability_, permeabilityLens_;
    Scalar eps_;
};

} // end namespace

#endif
