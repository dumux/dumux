// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Flow subproblem of the small-deformation Cryer benchmark: drained sphere surface,
 *        no flow across the symmetry planes.
 */
#ifndef DUMUX_TEST_POROELASTIC_CRYER_ONEP_PROBLEM_HH
#define DUMUX_TEST_POROELASTIC_CRYER_ONEP_PROBLEM_HH

#include <memory>
#include <string>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Flow subproblem of the small-deformation Cryer benchmark.
 */
template<class TypeTag>
class CryerOnePProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using GlobalPosition = typename GridDiscretization::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

public:
    CryerOnePProblem(std::shared_ptr<const GridDiscretization> gridDiscretization,
                     std::shared_ptr<GetPropType<TypeTag, Properties::SpatialParams>> spatialParams,
                     const std::string& paramGroup = "OneP")
    : ParentType(gridDiscretization, spatialParams, paramGroup)
    {
        problemName_ = getParam<std::string>("Vtk.OutputName") + "_"
            + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    const std::string& name() const
    { return problemName_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (isOnSymmetryPlane_(globalPos))
            values.setAllNeumann();
        else
            values.setAllDirichlet();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

    NumEqVector neumannAtPos(const GlobalPosition&) const
    { return NumEqVector(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

    NumEqVector sourceAtPos(const GlobalPosition&) const
    { return NumEqVector(0.0); }

private:
    bool isOnSymmetryPlane_(const GlobalPosition& globalPos) const
    {
        for (int d = 0; d < GlobalPosition::dimension; ++d)
            if (globalPos[d] < eps_)
                return true;
        return false;
    }

    static constexpr Scalar eps_ = 1e-6;
    std::string problemName_;
};

} // end namespace Dumux

#endif
