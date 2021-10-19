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
 * \ingroup Common
 * \brief Base class for all finite volume problems that are parameterized.
 */
#ifndef DUMUX_COMMON_FV_PROBLEM_WITH_SPATIAL_PARAMS_HH
#define DUMUX_COMMON_FV_PROBLEM_WITH_SPATIAL_PARAMS_HH

#include <memory>

#include <dumux/common/properties.hh>
#include <dumux/common/fvproblem.hh>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Base class for all finite-volume problems using spatial parameters.
 */
template<class TypeTag>
class FVProblemWithSpatialParams
: public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

public:
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    /*!
     * \brief Constructor
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     * \note This constructor assumes the spatial parameters to be constructible from a grid geometry
     */
    FVProblemWithSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                               const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , spatialParams_(std::make_shared<SpatialParams>(gridGeometry))
    {}

    /*!
     * \brief Constructor
     * \param gridGeometry The finite volume grid geometry
     * \param spatialParams The spatially varying parameters
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    FVProblemWithSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                               std::shared_ptr<SpatialParams> spatialParams,
                               const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , spatialParams_(spatialParams)
    {}

    //! Return a reference to the underlying spatial parameters
    const SpatialParams& spatialParams() const
    { return *spatialParams_; }

    //! Return a reference to the underlying spatial parameters
    SpatialParams& spatialParams()
    { return *spatialParams_; }

private:
    //! Spatially varying parameters
    std::shared_ptr<SpatialParams> spatialParams_;
};

} // end namespace Dumux

#endif
