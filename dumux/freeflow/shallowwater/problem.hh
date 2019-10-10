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
 * \ingroup ShallowWaterModel
 * \copydoc Dumux::ShallowWaterProblem
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_PROBLEM_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_PROBLEM_HH

#include <dumux/common/fvproblem.hh>
#include <dumux/common/properties.hh>
#include "model.hh"

namespace Dumux {

/*!
 * \ingroup ShallowWaterModel
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
