// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Mechanical subproblem of the small-deformation Cryer benchmark: pressure load on the
 *        sphere surface, zero normal displacement on the symmetry planes.
 */
#ifndef DUMUX_TEST_POROELASTIC_CRYER_POROELASTIC_PROBLEM_HH
#define DUMUX_TEST_POROELASTIC_CRYER_POROELASTIC_PROBLEM_HH

#include <algorithm>
#include <memory>
#include <string>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Mechanical subproblem of the small-deformation Cryer benchmark.
 *
 * The surface load is the traction \f$ \boldsymbol{\sigma} \mathbf{n} = -p_0 \mathbf{N} \f$ of the
 * total stress with the radial unit vector \f$ \mathbf{N} \f$.
 */
template<class TypeTag>
class CryerPoroElasticProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using GlobalPosition = typename GridDiscretization::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    static constexpr int dim = GlobalPosition::dimension;

public:
    CryerPoroElasticProblem(std::shared_ptr<const GridDiscretization> gridDiscretization,
                            std::shared_ptr<GetPropType<TypeTag, Properties::SpatialParams>> spatialParams,
                            const std::string& paramGroup = "PoroElastic")
    : ParentType(gridDiscretization, spatialParams, paramGroup)
    , p0_(getParam<Scalar>("Problem.PressureLoad"))
    {
        problemName_ = getParam<std::string>("Vtk.OutputName") + "_"
            + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    const std::string& name() const
    { return problemName_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        for (int d = 0; d < dim; ++d)
            if (globalPos[d] < eps_)
                values.setDirichlet(d);
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);
        if (std::none_of(globalPos.begin(), globalPos.end(), [](Scalar x){ return x < eps_; }))
        {
            values = globalPos;
            values *= -p0_/globalPos.two_norm();
        }
        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

private:
    static constexpr Scalar eps_ = 1e-6;
    Scalar p0_;
    std::string problemName_;
};

} // end namespace Dumux

#endif
