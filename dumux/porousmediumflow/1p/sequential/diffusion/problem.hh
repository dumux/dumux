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
 * \ingroup SequentialOnePModel
 * \brief Base class for all single phase diffusion problem
 */
#ifndef DUMUX_DIFFUSIONPROBLEM_1P_HH
#define DUMUX_DIFFUSIONPROBLEM_1P_HH

#include <dumux/porousmediumflow/sequential/onemodelproblem.hh>
#include <dumux/porousmediumflow/sequential/variableclass.hh>
#include <dumux/porousmediumflow/1p/sequential/diffusion/properties.hh>
#include <dumux/porousmediumflow/1p/sequential/celldata.hh>

namespace Dumux
{
namespace Properties
{
SET_TYPE_PROP(PressureOneP, Model, typename GET_PROP_TYPE(TypeTag, PressureModel));
}
/*!
 * \ingroup SequentialOnePModel
 * \brief  Base class for all single phase diffusion problem
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class DiffusionProblem1P: public OneModelProblem<TypeTag>
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, Problem);
    using ParentType = OneModelProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Grid = typename GridView::Grid;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    // material properties
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);

    using Element = typename GridView::Traits::template Codim<0>::Entity;

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief Constructs a DiffusionProblem1P object
     *
     * \param timeManager the time manager
     * \param gridView The grid view
     */
    DiffusionProblem1P(TimeManager& timeManager, Grid& grid)
    : ParentType(timeManager, grid), gravity_(0)
    {
        spatialParams_ = std::make_shared<SpatialParams>(asImp_());
        gravity_ = 0;
        if (getParam<bool>("Problem.EnableGravity"))
            gravity_[dim - 1] = -9.81;
    }

    /*!
     * \brief Constructs a DiffusionProblem1P object
     *
     * \param timeManager the time manager
     * \param gridView The grid view
     * \param spatialParams SpatialParams instantiation
     */
    DiffusionProblem1P(TimeManager& timeManager, Grid& grid, SpatialParams &spatialParams)
    : ParentType(timeManager, grid), gravity_(0)
    {
        spatialParams_ = Dune::stackobject_to_shared_ptr<SpatialParams>(spatialParams);
        gravity_ = 0;
        if (getParam<bool>("Problem.EnableGravity"))
            gravity_[dim - 1] = -9.81;
    }

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    DiffusionProblem1P(Grid& grid)
    : ParentType(grid, false), gravity_(0)
    {
        spatialParams_ = std::make_shared<SpatialParams>(asImp_());
        gravity_ = 0;
        if (getParam<bool>("Problem.EnableGravity"))
            gravity_[dim - 1] = -9.81;
    }

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     * \param spatialParams SpatialParams instantiation
     */
    DiffusionProblem1P(Grid& grid, SpatialParams& spatialParams)
    : ParentType(grid, false), gravity_(0)
    {
        spatialParams_ = Dune::stackobject_to_shared_ptr<SpatialParams>(spatialParams);
        gravity_ = 0;
        if (getParam<bool>("Problem.EnableGravity"))
            gravity_[dim - 1] = -9.81;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    //! \cond \private
    void timeIntegration()
    {
        //end simulation -> no time dependent problem!
        this->timeManager().setFinished();

        return;
    }

    void serialize()
    {
        return;
    }

    void deserialize(double t)
    {
        return;
    }
    //! \endcond

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param element The element
     *
     */
    Scalar temperature(const Element& element) const
    {
        return asImp_().temperatureAtPos(element.geometry().center());
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The position of the center of an element
     *
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
     *
     */
    Scalar referencePressure(const Element& element) const
    {
        return asImp_().referencePressureAtPos(element.geometry().center());
    }

    /*!
     * \brief Returns the reference pressure for evaluation of constitutive relations.
     *
     * \param globalPos The position of the center of an element
     *
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

    // \}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc SequentialOnePModel::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    GravityVector gravity_;

    // fluids and material properties
    std::shared_ptr<SpatialParams> spatialParams_;
};

}

#endif
