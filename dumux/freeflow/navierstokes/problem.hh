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

#warning "This header is deprecated and will be removed after release 3.6. Use the new mass and momemtum problem headers"

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/freeflow/navierstokes/staggered/problem.hh>

namespace Dumux {

//! The implementation is specialized for the different discretizations
template<class TypeTag, class DiscretizationMethod> struct NavierStokesParentProblemImpl;

// compatibility with old-style Navier-Stokes models
template<class TypeTag>
struct [[deprecated("Will be removed after 3.6. Directly use StaggeredFVProblem.")]]
NavierStokesParentProblemImpl<TypeTag, DiscretizationMethods::Staggered>
{
    using type = StaggeredFVProblem<TypeTag>;
};

//! The actual NavierStokesParentProblem
template<class TypeTag>
using NavierStokesParentProblem [[deprecated("Will be removed after 3.6. Directly use StaggeredFVProblem.")]]
= typename NavierStokesParentProblemImpl<
    TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod
>::type;

template<class TypeTag, class DiscretizationMethod>
class NavierStokesProblemImpl;

template<class TypeTag>
class [[deprecated("Class will be removed after 3.6. Use new staggered problem.")]]
NavierStokesProblemImpl<TypeTag, DiscretizationMethods::Staggered>
 : public NavierStokesStaggeredProblem<TypeTag>
{
    using ParentType = NavierStokesStaggeredProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

public:
    /*!
     * \brief The constructor for usage without a coupling manager
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    {}
};

template<class TypeTag>
class [[deprecated("Class will be removed after 3.6. Use matching momentum problem.")]]
NavierStokesProblemImpl<TypeTag, DiscretizationMethods::FCStaggered>
 : public NavierStokesMomentumProblemImpl<TypeTag, DiscretizationMethods::FCStaggered>
{
    using ParentType = NavierStokesMomentumProblemImpl<TypeTag, DiscretizationMethods::FCStaggered>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    //! These types are used in place of the typical NumEqVector type.
    //! In the numeqvector assembly type, only one equation per DOF (face) is considered
    //! while the type here provides one entry for each world dimension.
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;

    using MomentumFluxType = typename ParentType::MomentumFluxType;

    //! Export the boundary types.
    using BoundaryTypes = typename ParentType::BoundaryTypes;

    using ParentType::isMomentumProblem;

    /*!
     * \brief The constructor
     * \param gridGeometry The finite volume grid geometry
     * \param couplingManager The coupling manager (couples mass and momentum equations)
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            std::shared_ptr<CouplingManager> couplingManager,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, couplingManager, paramGroup)
    {}

    /*!
     * \brief The constructor for usage without a coupling manager
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    {}
};

template<class TypeTag>
class [[deprecated("Class will be removed after 3.6. Use matching momentum problem.")]]
NavierStokesProblemImpl<TypeTag, DiscretizationMethods::FCDiamond>
 : public NavierStokesMomentumProblemImpl<TypeTag, DiscretizationMethods::FCDiamond>
{
    using ParentType = NavierStokesMomentumProblemImpl<TypeTag, DiscretizationMethods::FCDiamond>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    //! These types are used in place of the typical NumEqVector type.
    //! In the numeqvector assembly type, only one equation per DOF (face) is considered
    //! while the type here provides one entry for each world dimension.
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;

    //! Export the boundary types.
    using BoundaryTypes = typename ParentType::BoundaryTypes;

    using ParentType::isMomentumProblem;

    /*!
     * \brief The constructor
     * \param gridGeometry The finite volume grid geometry
     * \param couplingManager The coupling manager (couples mass and momentum equations)
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            std::shared_ptr<CouplingManager> couplingManager,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, couplingManager, paramGroup)
    {}

    /*!
     * \brief The constructor for usage without a coupling manager
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    {}
};

template<class TypeTag>
class [[deprecated("Class will be removed after 3.6. Use matching mass problem.")]]
NavierStokesProblemImpl<TypeTag, DiscretizationMethods::CCTpfa>
 : public NavierStokesMassProblemImpl<TypeTag, DiscretizationMethods::CCTpfa>
{
    using ParentType = NavierStokesMassProblemImpl<TypeTag, DiscretizationMethods::CCTpfa>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    //! These types are used in place of the typical NumEqVector type.
    //! In the numeqvector assembly type, only one equation per DOF (face) is considered
    //! while the type here provides one entry for each world dimension.
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;

    //! Export the boundary types.
    using BoundaryTypes = typename ParentType::BoundaryTypes;

    using ParentType::isMomentumProblem;

    /*!
     * \brief The constructor
     * \param gridGeometry The finite volume grid geometry
     * \param couplingManager The coupling manager (couples mass and momentum equations)
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            std::shared_ptr<CouplingManager> couplingManager,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, couplingManager, paramGroup)
    {}

    /*!
     * \brief The constructor for usage without a coupling manager
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    {}
};

/*!
 * \ingroup NavierStokesModel
 * \brief Navier-Stokes problem class
 *
 * Inherit from this problem to implement Navier-Stokes problems
 */
template<class TypeTag>
using NavierStokesProblem [[deprecated("Will be removed after 3.6. Use matching mass or momentum problem.")]]
= NavierStokesProblemImpl<TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>;

} // end namespace Dumux

#endif
