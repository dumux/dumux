// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Base class for all finite volume problems that are parameterized.
 */
#ifndef DUMUX_COMMON_FV_PROBLEM_WITH_SPATIAL_PARAMS_HH
#define DUMUX_COMMON_FV_PROBLEM_WITH_SPATIAL_PARAMS_HH

#include <memory>

#include <dumux/common/properties.hh>
#include <dumux/common/fvproblem.hh>

namespace Dumux {

/*!
 * \ingroup Core
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
