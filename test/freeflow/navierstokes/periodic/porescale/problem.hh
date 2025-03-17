// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Pore Scale Simulations
 * \ingroup Basic Pore Structure
 * \brief The problem file for the evaluation of a pore geometry
 */

#ifndef DUMUX_PERIODIC_FLOW_PORE_TEST_PROBLEM_HH
#define DUMUX_PERIODIC_FLOW_PORE_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>

#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>

namespace Dumux {

template <class TypeTag, class BaseProblem>
class FlowPoreStructureProblem : public BaseProblem
{
    using ParentType = BaseProblem;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    FlowPoreStructureProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                             std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    , eps_(1e-6)
    {
        momentumForce_ = getParam<Scalar>("Problem.MomentumForce");
        referencePressure_ = 0.0;
        problemName_  =  "test_ff_periodic_subgrid";
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    { return problemName_ ; }

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    template<class ElementVolumeVariables>
    Sources source(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume& scv) const
    {
        Sources source(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            GlobalPosition globalPos = scv.dofPosition();
            if (isInlet_(globalPos))
            {
                const auto& frontalScvf = (*scvfs(fvGeometry, scv).begin());
                source[Indices::momentumXBalanceIdx] = momentumForce_ * frontalScvf.area() / scv.volume();
            }
        }

        return source;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    DirichletValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    { return DirichletValues(0.0); }

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is subtracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return referencePressure_; }

    // Return cellCenteredVelocities
    template <class SolutionVector>
    std::vector<VelocityVector> cellCenteredVelocity(const SolutionVector& curSol) const
    {
        std::vector<VelocityVector> velocity;
        velocity.resize(this->gridGeometry().elementMapper().size(), VelocityVector(0.0));

        // calculate cell-center-averaged velocities
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);

            // calculate velocities
            VelocityVector velocityTemp(0.0);
            for (auto&& scv : scvs(fvGeometry))
            {
                const int dofIdxFace = scv.dofIndex();
                const auto numericalSolutionFace = curSol[CouplingManager::freeFlowMomentumIndex]
                                                         [dofIdxFace]
                                                         [Indices::velocity(scv.dofAxis())];
                velocityTemp[scv.dofAxis()] += numericalSolutionFace;
            }

            for (unsigned int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                velocity[elementIdx][dimIdx] = velocityTemp[dimIdx] * 0.5; // faces are equidistant to cell center
        }
        return velocity;
    }

private:
    bool isInlet_(const GlobalPosition& globalPos) const
    { return (globalPos[0] <  this->gridGeometry().bBoxMin()[0] + eps_); }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_); }

    Scalar eps_;
    Scalar momentumForce_;
    Scalar referencePressure_;
    std::string problemName_;
};

} // end namespace Dumux

#endif
