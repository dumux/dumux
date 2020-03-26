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
 * \ingroup SequentialTwoPModel
 * \brief Base class for stationary solution of a two-phase diffusion/pressure equation.
 */
#ifndef DUMUX_DIFFUSIONPROBLEM_2P_HH
#define DUMUX_DIFFUSIONPROBLEM_2P_HH

#include <dumux/porousmediumflow/sequential/onemodelproblem.hh>
#include "properties.hh"

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief  Base class for stationary solution of a two-phase diffusion/pressure equation.
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class DiffusionProblem2P: public OneModelProblem<TypeTag>
{
    using Implementation = GetPropType<TypeTag, Properties::Problem>;
    using ParentType = OneModelProblem<TypeTag>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Grid = typename GridView::Grid;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;
    using PressureModel = GetPropType<TypeTag, Properties::PressureModel>;

    // material properties
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using Element = typename GridView::Traits::template Codim<0>::Entity;

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;

    // private!! copy constructor
    DiffusionProblem2P(const DiffusionProblem2P&)
    {}

public:
    /*!
     * \brief Constructs a DiffusionProblem2P object
     *
     * \param timeManager the time manager
     * \param grid The grid
     */
    DiffusionProblem2P(TimeManager& timeManager, Grid& grid)
    : ParentType(timeManager, grid), gravity_(0)
    {
        spatialParams_ = std::make_shared<SpatialParams>(asImp_());
        gravity_ = 0;
        if (getParam<bool>("Problem.EnableGravity"))
            gravity_[dim - 1] = -9.81;

        pressModel_ = std::make_shared<PressureModel>(asImp_());
    }
    /*!
     * \brief Constructs a DiffusionProblem2P object
     *
     * \param timeManager the time manager
     * \param grid The grid
     * \param spatialParams SpatialParams instantiation
     */
    DiffusionProblem2P(TimeManager& timeManager, Grid& grid, SpatialParams& spatialParams)
    : ParentType(timeManager, grid), gravity_(0)
    {
        spatialParams_ = Dune::stackobject_to_shared_ptr<SpatialParams>(spatialParams);
        gravity_ = 0;
        if (getParam<bool>("Problem.EnableGravity"))
            gravity_[dim - 1] = -9.81;

        pressModel_ = std::make_shared<PressureModel>(asImp_());
    }

    /*!
     * \brief Constructs a DiffusionProblem2P object
     *
     * \param grid The grid
     */
    DiffusionProblem2P(Grid& grid)
    : ParentType(grid, false), gravity_(0)
    {
        spatialParams_ = std::make_shared<SpatialParams>(asImp_());
        gravity_ = 0;
        if (getParam<bool>("Problem.EnableGravity"))
            gravity_[dim - 1] = -9.81;

        pressModel_ = std::make_shared<PressureModel>(asImp_());
    }
    /*!
     * \brief Constructs a DiffusionProblem2P object
     *
     * \param grid The grid
     * \param spatialParams SpatialParams instantiation
     */
    DiffusionProblem2P(Grid& grid, SpatialParams& spatialParams)
    : ParentType(grid, false), gravity_(0)
    {
        spatialParams_ = Dune::stackobject_to_shared_ptr<SpatialParams>(spatialParams);
        gravity_ = 0;
        if (getParam<bool>("Problem.EnableGravity"))
            gravity_[dim - 1] = -9.81;

        pressModel_ = std::make_shared<PressureModel>(asImp_());
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Time integration function called by the time manager
     *
     * For stationary diffusion problems this function just finishes the simulation.
     */
    void timeIntegration()
    {
        //end simulation -> no time dependent problem!
        this->timeManager().setFinished();

        return;
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param element The element
     */
    Scalar temperature(const Element& element) const
    {
        return this->asImp_().temperatureAtPos(element.geometry().center());
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The position of the center of an element
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        // Throw an exception (there is no initial condition)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a temperatureAtPos() method.");
    }

    /*!
     * \brief Returns the reference pressure for evaluation of constitutive relations.
     *
     * \param element The element
     */
    Scalar referencePressure(const Element& element) const
    {
        return this->asImp_().referencePressureAtPos(element.geometry().center());
    }

    /*!
     * \brief Returns the reference pressure for evaluation of constitutive relations.
     *
     * \param globalPos The position of the center of an element
     */
    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
        // Throw an exception (there is no initial condition)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a referencePressureAtPos() method.");
    }

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>EnableGravity</tt> property is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GravityVector &gravity() const
    {
        return gravity_;
    }

    /*!
     * \brief Returns the spatial parameters object.
     */
    SpatialParams &spatialParams()
    {
        return *spatialParams_;
    }

    /*!
     * \brief Returns the spatial parameters object.
     */
    const SpatialParams &spatialParams() const
    {
        return *spatialParams_;
    }

    /*!
     * \brief Returns the pressure model used for the problem.
     */
    PressureModel &pressureModel()
    { return *pressModel_; }

    //! \copydoc IMPETProblem::pressureModel() // TODO doc me!
    const PressureModel &pressureModel() const
    { return *pressModel_; }

    // \}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc IMPETProblem::asImp_() // TODO doc me!
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    GravityVector gravity_;

    // fluids and material properties
    std::shared_ptr<SpatialParams> spatialParams_;
    bool newSpatialParams_;
    std::shared_ptr<PressureModel> pressModel_;
};

} // end namespace Dumux

#endif
