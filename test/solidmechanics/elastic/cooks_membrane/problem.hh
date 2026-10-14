// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief Cook's membrane: a tapered panel clamped on the left edge and loaded by a shear traction
 *        on the right edge (plane strain linear elasticity).
 *
 * The panel has the corners (0, 0), (48, 44), (48, 60) and (0, 44). The left edge \f$ x = 0 \f$ is
 * clamped, the right edge \f$ x = 48 \f$ carries the traction \f$ (0, t) \f$, and the other edges are
 * traction-free.
 */
#ifndef DUMUX_TEST_ELASTIC_COOKS_MEMBRANE_PROBLEM_HH
#define DUMUX_TEST_ELASTIC_COOKS_MEMBRANE_PROBLEM_HH

#include <vector>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

template<class TypeTag>
class CooksMembraneProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<PrimaryVariables::size()>;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, PrimaryVariables,
                                                                   typename IndexTraits<GridView>::GridIndex>;

public:
    CooksMembraneProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , shearTraction_(getParam<Scalar>("Problem.ShearTraction"))
    { appendDirichletConstraints_(); }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    //! The clamped left edge is constrained, all other edges carry a (possibly zero) traction
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] > eps_)
            values.setAllFluxBoundary();
        return values;
    }

    //! The traction on the right edge, \f$ \mathbf{g} = -\boldsymbol{\sigma}\mathbf{n} = (0, -t) \f$
    NumEqVector boundaryFluxAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector flux(0.0);
        if (globalPos[0] > 48.0 - eps_)
            flux[1] = -shearTraction_;
        return flux;
    }

private:
    //! The local dofs of the clamped edge have zero displacement
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
                    constraints_.push_back(
                        DirichletConstraintData{std::move(info), PrimaryVariables(0.0), localDof.dofIndex()}
                    );
                }
            }
        }
    }

    static constexpr Scalar eps_ = 1e-6;
    Scalar shearTraction_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
