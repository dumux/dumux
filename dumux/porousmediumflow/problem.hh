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
 * \ingroup PorousmediumflowModels
 * \brief Base class for all porous media problems.
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_PROBLEM_HH
#define DUMUX_POROUS_MEDIUM_FLOW_PROBLEM_HH

#include <dumux/common/fvproblem.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief Base class for all fully implicit porous media problems.
 *
 * TODO: derive from base problem property?
 */
template<class TypeTag>
class PorousMediumFlowProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;


public:
    //! Export spatial parameter type
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    /*!
     * \brief Constructor, passing the spatial parameters.
     *
     * \param gridGeometry The finite volume grid geometry
     * \param spatialParams The spatial parameter class
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    PorousMediumFlowProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                            std::shared_ptr<SpatialParams> spatialParams,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , gravity_(0.0)
    , spatialParams_(spatialParams)
    {
        const bool enableGravity = getParamFromGroup<bool>(paramGroup, "Problem.EnableGravity");
        if (enableGravity)
            gravity_[dimWorld-1]  = -9.81;
    }

    /*!
     * \brief Constructor, constructing the spatial parameters.
     *
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    PorousMediumFlowProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                            const std::string& paramGroup = "")
    : PorousMediumFlowProblem(gridGeometry,
                              std::make_shared<SpatialParams>(gridGeometry),
                              paramGroup)
    {}

    /*!
     * \name Physical parameters for porous media problems
     */
    // \{

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at a given global position.
     *
     * This is not specific to the discretization. By default it just
     * calls temperature().
     *
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    { return this->asImp_().temperature(); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    {
        DUNE_THROW(Dune::NotImplemented, "temperature() method not implemented by the user problem");
    }

    /*!
     * \brief Returns the spatial parameters object.
     */
    SpatialParams &spatialParams()
    { return *spatialParams_; }

    /*!
     * \brief Returns the spatial parameters object.
     */
    const SpatialParams &spatialParams() const
    { return *spatialParams_; }

    // \}

protected:
    //! The gravity acceleration vector
    GravityVector gravity_;

    // material properties of the porous medium
    std::shared_ptr<SpatialParams> spatialParams_;
};

} // end namespace Dumux

#endif
