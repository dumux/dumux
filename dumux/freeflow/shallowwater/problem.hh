// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterModels
 * \copydoc Dumux::ShallowWaterProblem
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_PROBLEM_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_PROBLEM_HH

#include <dumux/common/fvproblem.hh>
#include <dumux/common/properties.hh>
#include "model.hh"

namespace Dumux {

/*!
 * \ingroup ShallowWaterModels
 * \brief Shallow water problem base class.
 */
template<class TypeTag>
class ShallowWaterProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

public:
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    /*!
     * \brief Constructor, passing the spatial parameters
     *
     * \param gridGeometry The finite volume grid geometry
     * \param spatialParams The spatial parameter class
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    ShallowWaterProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                        std::shared_ptr<SpatialParams> spatialParams,
                        const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , spatialParams_(spatialParams)
    {}

    /*!
     * \brief Constructor, constructing the spatial parameters
     *
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    ShallowWaterProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                        const std::string& paramGroup = "")
    : ShallowWaterProblem(gridGeometry,
                          std::make_shared<SpatialParams>(gridGeometry),
                          paramGroup)
    {}



    /*!
     * \brief Returns the spatial parameters object.
     */
    const SpatialParams &spatialParams() const
    { return *spatialParams_; }

    // \}


private:
    std::shared_ptr<SpatialParams> spatialParams_; //!< the spatial parameters
};

} // end namespace Dumux

#endif
