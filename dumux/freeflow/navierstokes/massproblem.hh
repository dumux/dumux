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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesProblem
 */
#ifndef DUMUX_NAVIERSTOKES_PROBLEM_HH
#define DUMUX_NAVIERSTOKES_PROBLEM_HH

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/momentum/boundarytypes.hh>

namespace Dumux {

//! The implementation is specialized for the different discretizations
template<class TypeTag, class DiscretizationMethod> struct NavierStokesParentProblemImpl;

//! The actual NavierStokesParentProblem
template<class TypeTag>
using NavierStokesParentProblem = typename NavierStokesParentProblemImpl<
    TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod
>::type;

template<class TypeTag, class DiscretizationMethod>
class NavierStokesProblemImpl;

template<class TypeTag>
class NavierStokesProblemImpl<TypeTag, DiscretizationMethods::CCTpfa>
: public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;
    using VelocityVector = GlobalPosition;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    //using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr bool isCoupled_ = false;//!std::is_empty_v<CouplingManager>;

public:

    //! These types are used in place of the typical NumEqVector type.
    //! In the momentum problem, there is a descrepancy between
    //! the assembly NumEqVector type, and the type used in the residual.
    //! These aliases are used in both problems to differentiate.
    using InitialValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using Sources = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using DirichletValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using BoundaryFluxes = Dune::FieldVector<Scalar, ModelTraits::numEq()>;

    //! Export the boundary types.
    using BoundaryTypes = Dumux::BoundaryTypes<ModelTraits::numEq()>;

    //! this problem is used for the mass balance model
    static constexpr bool isMomentumProblem() { return false; }

    /*!
     * \brief The constructor
     * \param gridGeometry The finite volume grid geometry
     * \param couplingManager The coupling manager (couples mass and momentum equations)
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    {}

    /*!
     * \brief The constructor for usage without a coupling manager
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            const std::string& paramGroup = "")
    : NavierStokesProblemImpl(gridGeometry, {}, paramGroup)
    {}

    /*!
     * \brief Returns the normal velocity at a given sub control volume face.
     */
    VelocityVector faceVelocity(const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const SubControlVolumeFace& scvf) const
    {
        if constexpr (isCoupled_)
            return couplingManager_->faceVelocity(element, scvf);
        else
            return asImp_().velocityAtPos(scvf.ipGlobal());
    }

    /*!
     * \brief Returns the velocity at a given position.
     */
    VelocityVector velocityAtPos(const GlobalPosition&) const
    {
        DUNE_THROW(Dune::NotImplemented, "velocityAtPos not implemented");
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at a given global position.
     *
     * This is not specific to the discretization. By default it just
     * calls temperature().
     *
     * \param globalPos The position in global coordinates where the temperature should be specified.
     * \param deprecationHelper The deprecation helper
     */
    [[deprecated("temperature should now be defined in the spatial params with temperature(globalPos)")]]
    Scalar temperatureAtPos(const GlobalPosition& globalPos, int deprecationHelper = 0) const
    { return asImp_().temperature(); }

    /*!
     * \brief Returns the temperature within the domain.
     */
    [[deprecated("temperature should now be defined in the spatial params with temperature(globalPos)")]]
    Scalar temperature(int deprecationHelper = 0) const
    {
        DUNE_THROW(Dune::NotImplemented, "temperature() method not implemented by the actual problem");
    }

    const CouplingManager& couplingManager() const
    {
        if constexpr (isCoupled_)
            return *couplingManager_;
        else
            DUNE_THROW(Dune::InvalidStateException,
                "Accessing coupling manager of an uncoupled problem is not possible."
            );
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    std::shared_ptr<CouplingManager> couplingManager_;
};

/*!
 * \ingroup NavierStokesModel
 * \brief Navier-Stokes problem class
 *
 * Inherit from this problem to implement Navier-Stokes problems
 */
template<class TypeTag>
using NavierStokesProblem = NavierStokesProblemImpl<
    TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod
>;

} // end namespace Dumux

#endif
