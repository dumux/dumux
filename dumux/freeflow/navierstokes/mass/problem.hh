// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMassProblem
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_PROBLEM_HH
#define DUMUX_NAVIERSTOKES_MASS_PROBLEM_HH

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

// default implementation
template<class TypeTag, class DiscretizationMethod>
class NavierStokesMassProblemImpl : public FVProblemWithSpatialParams<TypeTag>
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
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr bool isCoupled_ = !std::is_empty_v<CouplingManager>;

public:

    //! These types are used in place of the typical NumEqVector type.
    //! In the momentum problem, there is a discrepancy between
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
    NavierStokesMassProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                                std::shared_ptr<CouplingManager> couplingManager,
                                const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , couplingManager_(couplingManager)
    {}

    /*!
     * \brief The constructor for usage without a coupling manager
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesMassProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                                const std::string& paramGroup = "")
    : NavierStokesMassProblemImpl(gridGeometry, {}, paramGroup)
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
     * \brief Returns the velocity at the element center.
     */
    VelocityVector elementVelocity(const FVElementGeometry& fvGeometry) const
    {
        if constexpr (isCoupled_)
            return couplingManager_->elementVelocity(fvGeometry);
        else
            return asImp_().velocityAtPos(fvGeometry.element().geometry().center());
    }

    /*!
     * \brief Returns the velocity at a given position.
     */
    VelocityVector velocityAtPos(const GlobalPosition&) const
    {
        DUNE_THROW(Dune::NotImplemented, "velocityAtPos not implemented");
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
 * \brief Navier-Stokes mass problem class
 *
 * Inherit from this problem to implement Navier-Stokes mass problems
 */
template<class TypeTag>
using NavierStokesMassProblem = NavierStokesMassProblemImpl<
    TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod
>;

} // end namespace Dumux

#endif
