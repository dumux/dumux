// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Fluid pressure subdomain problem for the Cryer sphere benchmark.
 *
 * Boundary conditions:
 * - Sphere surface (drained): p_f = 0 (constraint)
 * - Symmetry planes: no flow
 */
#ifndef DUMUX_CRYER_PROBLEM_FLUID_PRESSURE_HH
#define DUMUX_CRYER_PROBLEM_FLUID_PRESSURE_HH

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

template<class TypeTag>
class CryerFluidPressureProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using BoundaryFace = typename GridGeometry::BoundaryFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using ConstraintValues = Dune::FieldVector<Scalar, numEq>;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

public:
    CryerFluidPressureProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                               std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    { appendDirichletConstraints_(); }

    const CouplingManager& couplingManager() const { return *couplingManager_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& pos) const
    {
        BoundaryTypes bct;
        if (isOnSymmetryPlane_(pos))
            bct.setAllFluxBoundary();
        return bct;
    }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

    const auto& constraints() const { return constraints_; }

private:
    //! The fluid pressure is zero on the drained sphere surface
    void appendDirichletConstraints_()
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                if (this->boundaryTypes(elemDisc, boundaryFace).hasOnlyFluxBoundary())
                    continue;

                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    ConstraintInfo info;
                    info.setAll();
                    constraints_.push_back(DirichletConstraintData{std::move(info), ConstraintValues(0.0), localDof.dofIndex()});
                }
            }
        }
    }

    bool isOnSymmetryPlane_(const GlobalPosition& pos) const
    {
        for (int d = 0; d < GridView::dimension; ++d)
            if (std::abs(pos[d]) < eps_) return true;
        return false;
    }

    static constexpr Scalar eps_ = 1e-8;
    std::shared_ptr<CouplingManager> couplingManager_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux
#endif
