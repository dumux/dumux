// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Base class for all problems using spatial parameters.
 */
#ifndef DUMUX_COMMON_PROBLEM_WITH_SPATIAL_PARAMS_HH
#define DUMUX_COMMON_PROBLEM_WITH_SPATIAL_PARAMS_HH

#include <memory>

#include <dumux/common/properties.hh>
#include <dumux/common/problem.hh>

namespace Dumux {

/*!
 * \ingroup Core
 * \brief Base class for all problems using spatial parameters.
 */
template<class TypeTag, class BaseProblem = Problem<TypeTag>>
class ProblemWithSpatialParams
: public BaseProblem
{
    using ParentType = BaseProblem;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

public:
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    /*!
     * \brief Constructor
     * \param gridGeometry The discretization related grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     * \note This constructor assumes the spatial parameters to be constructible from a grid geometry
     */
    ProblemWithSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                             const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , spatialParams_(std::make_shared<SpatialParams>(gridGeometry))
    {}

    /*!
     * \brief Constructor
     * \param gridGeometry The discretization related grid geometry
     * \param spatialParams The spatially varying parameters
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    ProblemWithSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
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
