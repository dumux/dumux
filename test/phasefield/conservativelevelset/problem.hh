// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ConservativeLevelSetTests
 * \brief Test problem for the conservative level-set model discretized with
 *        PQ2 (hybrid CVFE) elements.
 */
#ifndef DUMUX_TEST_PHASEFIELD_CONSERVATIVE_LEVEL_SET_PROBLEM_HH
#define DUMUX_TEST_PHASEFIELD_CONSERVATIVE_LEVEL_SET_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problem.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux {

/*!
 * \ingroup ConservativeLevelSetTests
 * \brief Test problem for the conservative level-set model on PQ2 (hybrid
 *        CVFE) elements.
 *
 * The problem sets homogeneous Neumann (no-flux) boundary conditions on the
 * entire domain boundary and has no source term; the reinitialization
 * equation only redistributes phi to keep a fixed-width interface profile.
 */
template<class TypeTag>
class ConservativeLevelSetTestProblem : public Dumux::Experimental::Problem<TypeTag>
{
    using ParentType = Dumux::Experimental::Problem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<ModelTraits::numEq()>;

public:
    ConservativeLevelSetTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        epsilon_ = getParam<Scalar>("Problem.InterfaceThickness");
    }

    /*!
     * \brief Specifies which kind of boundary condition should be used for
     *        which equation on a given boundary position. We use a
     *        homogeneous Neumann (no-flux) condition everywhere.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();
        return values;
    }

    //! the target interface half-width of the reinitialized profile
    Scalar epsilon() const
    { return epsilon_; }

private:
    Scalar epsilon_;
};

} // end namespace Dumux

#endif
