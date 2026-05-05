// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_MEMBRANE_PLATE_TEST_PROBLEM_HH
#define DUMUX_MEMBRANE_PLATE_TEST_PROBLEM_HH

#include <cmath>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup MembranePlateModel
 * \brief Test problem: clamped circular membrane under uniform load
 *
 * Solves \f$ -T\,\nabla^2 w = F \f$ on a disk of radius \f$ R \f$
 * with homogeneous Dirichlet conditions \f$ w = 0 \f$ on \f$ \partial\Omega \f$.
 * The analytic solution is
 * \f[ w(r) = \frac{F}{4T}(R^2 - r^2). \f]
 */
template<class TypeTag>
class MembranePlateTestProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

public:
    MembranePlateTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , force_(getParam<Scalar>("Problem.Force"))
    , tension_(getParam<Scalar>("Problem.Tension"))
    , radius_(getParam<Scalar>("Problem.Radius"))
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector source(0.0);
        source[Indices::deformationEqIdx] = force_;
        return source;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! Analytic solution: clamped circular membrane under uniform load
    Scalar analyticDeformation(const GlobalPosition& globalPos) const
    {
        const auto r2 = globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1];
        return force_/(4*tension_) * (radius_*radius_ - r2);
    }

private:
    Scalar force_, tension_, radius_;
};

} // end namespace Dumux

#endif
