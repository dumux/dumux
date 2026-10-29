// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Single-phase flow through two layers of different permeability
 */
#ifndef DUMUX_ONEP_LAYERED_TEST_PROBLEM_HH
#define DUMUX_ONEP_LAYERED_TEST_PROBLEM_HH

#include <memory>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Single-phase flow through two layers of different permeability
 *
 * The flow is in the x-direction with unit Darcy velocity in both layers. The pressure is
 * linear in each layer, and its gradient jumps at the interface by the ratio of the permeabilities.
 * The whole boundary carries the exact pressure.
 */
template<class TypeTag>
class OnePLayeredProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridDiscretization::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<ModelTraits::numEq()>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using ConstraintValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

public:
    OnePLayeredProblem(std::shared_ptr<const GridDiscretization> gridDiscretization)
    : ParentType(gridDiscretization)
    { appendDirichletConstraints_(); }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    { return BoundaryTypes{}; }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    //! The exact pressure, zero at the interface
    Scalar exactPressure(const GlobalPosition& globalPos) const
    {
        const auto& spatialParams = this->spatialParams();
        const Scalar distance = globalPos[0] - spatialParams.interfacePosition();
        return distance < 0.0 ? -distance/spatialParams.permeabilityLeft() : -distance/spatialParams.permeabilityRight();
    }

    //! The Dirichlet values the assembler constrains the solution to
    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

private:
    void appendDirichletConstraints_()
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    ConstraintInfo info;
                    info.setAll();
                    ConstraintValues values(exactPressure(ipData(elemDisc, localDof).global()));
                    constraints_.push_back(
                        DirichletConstraintData{std::move(info), std::move(values), localDof.dofIndex()}
                    );
                }
        }
    }

    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
