// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The problem for the DFG benchmark of the laminar flow around a cylinder (Navier-Stokes)
 *        using the new CVFE problem and variables interfaces
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_FLOW_AROUND_CYLINDER_PROBLEM_NEWINTERFACE_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_FLOW_AROUND_CYLINDER_PROBLEM_NEWINTERFACE_HH

#include <array>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/typetraits/griddiscretization.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>

namespace Dumux {

/*!
 * \brief Stationary laminar flow around a cylinder in a 2D channel (case 2D-1, Re = 20)
 *        of the DFG benchmark \cite SchaeferTurek1996 using the new CVFE problem interface
 *
 * The boundary conditions are the same as in problem.hh (old problem interface). The Dirichlet
 * conditions are imposed as constraints, and the drag and lift forces are computed by integrating
 * the traction over the cylinder faces with a quadrature rule of order 5.
 */
template <class TypeTag, class BaseProblem>
class DFGChannelTestProblemNewInterface : public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using DirichletValues = typename ParentType::DirichletValues;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using ConstraintValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<typename GridDiscretization::GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

    static constexpr int dim = GridDiscretization::GridView::dimension;
    static constexpr auto dimWorld = GridDiscretization::GridView::dimensionworld;

    using Element = typename GridDiscretization::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    DFGChannelTestProblemNewInterface(std::shared_ptr<const GridDiscretization> gridDiscretization,
                                      std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridDiscretization, couplingManager, ParentType::isMomentumProblem() ? "Momentum" : "Mass")
    {
        density_ = getParam<Scalar>("Component.LiquidDensity");
        viscosity_ = getParam<Scalar>("Component.LiquidDynamicViscosity");

        for (int i = 0; i < dimWorld; ++i)
            domainSize_[i] = gridDiscretization->bBoxMax()[i] - gridDiscretization->bBoxMin()[i];

        updateConstraints_();
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on the boundary.
     *
     * \param globalPos The position at the boundary
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            if (isOutlet_(globalPos))
                values.setAllFluxBoundary();
        }
        else
            values.setFluxBoundary(Indices::conti0EqIdx);

        return values;
    }

    /*!
     * \brief Evaluates the Dirichlet values (inflow profile and no-slip) at a given position
     */
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    {
        if (isInlet_(globalPos))
            return inflowVelocity_(globalPos[1]);
        else
            return DirichletValues(0.0);
    }

    /*!
     * \brief Return the Dirichlet boundary constraints
     */
    const auto& constraints() const
    { return constraints_; }

    /*!
     * \brief Evaluates the boundary flux related to a localDof at a given interpolation point.
     *
     * \param elemDisc The element discretization
     * \param elemVars All variables related to the element
     * \param faceIpData Interpolation point data
     */
    template<class ElementVariables, class FaceIpData>
    BoundaryFluxes boundaryFlux(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const FaceIpData& faceIpData) const
    {
        BoundaryFluxes values(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            if (this->enableInertiaTerms() && isOutlet_(faceIpData.global()))
            {
                // advective term: vv*n
                const auto& element = elemDisc.element();
                const auto elemSol = elementSolution(element, elemVars, elemDisc);
                const auto v = evalSolutionAtLocalPos(element, elemDisc.elementGeometry(), elemDisc.gridDiscretization(), elemSol, faceIpData.local());
                values.axpy(density_*(v*faceIpData.unitOuterNormal()), v);
            }
        }
        else
        {
            const auto& scvf = elemDisc.scvf(faceIpData.scvfIndex());
            values[Indices::conti0EqIdx] = this->velocity(elemDisc, faceIpData) * density_ * scvf.unitOuterNormal();
        }

        return values;
    }

    template<class ScvOrScvfOrIpData>
    Scalar density(const Element&,
                   const ElementDiscretization&,
                   const ScvOrScvfOrIpData&,
                   const bool isPreviousTimeStep = false) const
    {
        return density_;
    }

    template<class ScvOrScvfOrIpData>
    Scalar effectiveViscosity(const Element&,
                              const ElementDiscretization&,
                              const ScvOrScvfOrIpData&,
                              const bool isPreviousTimeStep = false) const
    {
        return viscosity_;
    }

    //! Computes pressure difference benchmark indicator
    //! (only works for correct domain size 2.2 x 0.41 and boundary conditions, Re=20)
    template<class SolutionVector, class GridVariables, class P = ParentType, typename std::enable_if_t<!P::isMomentumProblem(), int> = 0>
    Scalar evalPressureDifference(const GridVariables&, const SolutionVector& p) const
    {
        const auto& gridDiscretization = Dumux::gridDiscretization(*this);
        const auto& tree = gridDiscretization.boundingBoxTree();
        GlobalPosition evalPoint1({0.15, 0.2});
        GlobalPosition evalPoint2({0.25, 0.2});
        const Scalar stepSize = 0.001;

        // the evaluation points lie on the cylinder surface, so we may have to step into the fluid domain
        const auto evalPressure = [&](auto pos, bool backwards)
        {
            auto entities = intersectingEntities(pos, tree);
            while (entities.empty())
            {
                pos[0] += backwards ? -stepSize : stepSize;
                entities = intersectingEntities(pos, tree);
            }

            const auto element = gridDiscretization.element(entities[0]);
            const auto elemSol = elementSolution(element, p, gridDiscretization);
            return evalSolution(element, element.geometry(), gridDiscretization, elemSol, pos)[Indices::pressureIdx];
        };

        return evalPressure(evalPoint1, true) - evalPressure(evalPoint2, false);
    }

    //! Computes drag and lift coefficient benchmark indicator
    //! (only works for correct domain size 2.2 x 0.41 and boundary conditions, Re=20)
    template<class SolutionVector, class GridVariables, class P = ParentType, typename std::enable_if_t<P::isMomentumProblem(), int> = 0>
    auto evalDragAndLiftCoefficient(const GridVariables&, const SolutionVector& v) const
    {
        const auto& gridDiscretization = Dumux::gridDiscretization(*this);
        auto elemDisc = localView(gridDiscretization);

        // integrate the traction (mu grad(v) - pI)n over the cylinder surface
        BoundaryFluxes forces(0.0);
        Scalar cylinderSurface = 0.0;
        for (const auto& element : elements(gridDiscretization.gridView()))
        {
            elemDisc.bind(element);
            if (!elemDisc.hasBoundaryFaces())
                continue;

            const auto elemSol = elementSolution(element, v, gridDiscretization);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                if (!onCylinderBoundary_(boundaryFace.center()))
                    continue;

                using namespace Dumux::QuadratureRules;
                for (const auto& qpData : CVFE::quadratureRule(elemDisc, boundaryFace, DuneQuadrature<5>{}))
                {
                    const auto& ipData = qpData.ipData();
                    const auto& normal = ipData.unitOuterNormal();
                    const auto gradV = evalGradientsAtLocalPos(element, elemDisc.elementGeometry(), gridDiscretization, elemSol, ipData.local());
                    const auto pressure = this->pressure(element, elemDisc, ipData);

                    BoundaryFluxes traction(0.0);
                    for (int dir = 0; dir < dim; ++dir)
                        traction[dir] = viscosity_*(gradV[dir]*normal) - pressure*normal[dir];

                    forces.axpy(qpData.weight(), traction);
                    cylinderSurface += qpData.weight();
                }
            }
        }

        const Scalar dragForce = -forces[0];
        const Scalar liftForce = -forces[1];

        const Scalar uMean = 0.2;
        const Scalar lChar = 0.1;
        const Scalar coefficientFactor = 2.0/(uMean*uMean*lChar);

        std::cout << "Reynolds number: " << uMean*lChar/viscosity_ << std::endl;
        std::cout << "Cylinder surface: " << cylinderSurface
                  << " (reference: 0.31415926535)"
                  << std::endl;

        return std::make_pair( coefficientFactor*dragForce, coefficientFactor*liftForce );
    }

private:
    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridDiscretization().bBoxMin()[0] + eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridDiscretization().bBoxMax()[0] - eps_; }

    bool onCylinderBoundary_(const GlobalPosition& globalPos) const
    { return std::hypot(globalPos[0] - 0.2, globalPos[1] - 0.2) < 0.06; }

    DirichletValues inflowVelocity_(Scalar y) const
    {
        constexpr Scalar maxVelocity = 0.3;
        return { 4*maxVelocity*y*(domainSize_[1]-y)/(domainSize_[1]*domainSize_[1]), 0.0 };
    }

    void updateConstraints_() const
    {
        constraints_.clear();
        if constexpr (ParentType::isMomentumProblem())
            appendDirichletConstraints_();
    }

    void appendDirichletConstraints_() const
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                if (isOutlet_(boundaryFace.center()))
                    continue;

                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    // the corner dofs at the outlet belong to the outlet (do-nothing condition)
                    const auto dofPosition = ipData(elemDisc, localDof).global();
                    if (isOutlet_(dofPosition))
                        continue;

                    ConstraintInfo info;
                    info.setAll();
                    ConstraintValues values(dirichletAtPos(dofPosition));
                    constraints_.push_back(DirichletConstraintData{std::move(info), std::move(values), localDof.dofIndex()});
                }
            }
        }
    }

    static constexpr Scalar eps_ = 1e-10;
    Scalar density_, viscosity_;
    std::array<Scalar, dimWorld> domainSize_;
    mutable std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
