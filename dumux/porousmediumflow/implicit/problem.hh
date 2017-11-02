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
 * \brief Base class for all fully implicit porous media problems
 */
#ifndef DUMUX_IMPLICIT_POROUS_MEDIA_PROBLEM_HH
#define DUMUX_IMPLICIT_POROUS_MEDIA_PROBLEM_HH

#include <dumux/implicit/problem.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitBaseProblems
 * \brief Base class for all fully implicit porous media problems
 */
template<class TypeTag>
class ImplicitPorousMediaProblem : public ImplicitProblem<TypeTag>
{
    using ParentType = ImplicitProblem<TypeTag>;

    using Implementation = typename GET_PROP_TYPE(TypeTag, Problem);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::Traits::template Codim<dim>::Entity;

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     * \param verbose Turn verbosity on or off
     */
    ImplicitPorousMediaProblem(TimeManager &timeManager,
                const GridView &gridView,
                const bool verbose = true)
        : ParentType(timeManager, gridView),
          gravity_(0)
    {
        spatialParams_ = std::make_shared<SpatialParams>(asImp_(), gridView);

        if (getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.AddVelocity"))
            gravity_[dimWorld-1]  = -9.81;
    }

    /*!
     * \brief Called by the TimeManager in order to
     *        initialize the problem.
     *
     * We initialize the spatial parameters here.
     */
    void init()
    {
        ParentType::init();
        spatialParams_->init();
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Evaluate the initial phase state at an element (for cc models) or vertex (for box models)
     *
     * \param entity The entity (element or vertex)
     */
    template<class Entity>
    int initialPhasePresence(const Entity& entity)
    {
        static_assert(int(Entity::codimension) == 0 || int(Entity::codimension) == dim, "Entity must be element or vertex");

        return asImp_().initialPhasePresenceAtPos(entity.geometry().center());
    }

    /*!
     * \brief Evaluate the initial phase state at a given position
     *
     * \param globalPos The global position
     */
    int initialPhasePresenceAtPos(const GlobalPosition &globalPos) const
    {
       DUNE_THROW(Dune::InvalidStateException,
                  "The problem does not provide a initialPhasePresence() "
                   "or initialPhasePresenceAtPos() method.");
       return 0;
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at a given global position.
     *
     * This is not specific to the discretization. By default it just
     * calls temperature().
     *
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    { return asImp_().temperature(); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    { DUNE_THROW(Dune::NotImplemented, "temperature() method not implemented by the actual problem"); }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This is discretization independent interface. By default it
     * just calls gravity().
     */
    const GlobalPosition &gravityAtPos(const GlobalPosition &pos) const
    { return asImp_().gravity(); }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This method is used for problems where the gravitational
     * acceleration does not depend on the spatial position. The
     * default behaviour is that if the <tt>ProblemEnableGravity</tt>
     * property is true, \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$ holds,
     * else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$.
     */
    const GlobalPosition &gravity() const
    { return gravity_; }

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
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    GlobalPosition gravity_;

    // fluids and material properties
    std::shared_ptr<SpatialParams> spatialParams_;
};

}

#endif
