// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \ingroup TwoPTests
 * \brief The properties for the incompressible 2p test.
 */
#ifndef DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROBLEM_HH
#define DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/droplet/dropsolver.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief The incompressible 2p test problem.
 */
template<class TypeTag>
class TwoPTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using DropSolver = DropletSolverTwoP<TypeTag, false>;
    enum {
        waterPressureIdx = Indices::pressureIdx,
        airSaturationIdx = Indices::saturationIdx,
        contiDNAPLEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx,
        waterPhaseIdx = FluidSystem::phase0Idx,
        airPhaseIdx = FluidSystem::phase1Idx
    };

public:
    TwoPTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, GlobalPosition tabletCenter, Scalar tabletRadius)
    : ParentType(gridGeometry)
    , tabletCenter_(tabletCenter)
    , tabletRadius_{tabletRadius}
    {}

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (onCircleBoundary(globalPos))
        {
            bcTypes.setAllDirichlet();
        }
        else if (onUpperBoundary(globalPos))
        {
            if (dropletSolver_->isCoupledWithDroplet(globalPos))
            {
                bcTypes.setAllDirichlet();
            }
            else
            {
                bcTypes.setAllNeumann();
            }
        }
        else
            bcTypes.setAllNeumann();
        return bcTypes;
    }

    // auto boundaryTypes(const Element &element,
    //                    const SubControlVolume &scv) const
    // {
    //     BoundaryTypes bcTypes;
    //     const auto& globalPos = scv.dofPosition();

    //     if (onLeftBoundary(globalPos) || onRightBoundary(globalPos))
    //     {
    //         bcTypes.setAllDirichlet();
    //     }
    //     else if (onUpperBoundary(globalPos))
    //     {
    //         if (dropletSolver_->isCoupledWithDroplet(globalPos))
    //         {
    //             bcTypes.setAllDirichlet();
    //         }
    //         else
    //         {
    //             bcTypes.setAllNeumann();
    //         }
    //     }
    //     else
    //         bcTypes.setAllNeumann();
    //     return bcTypes;
    // }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The global position
     */
    // PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    // {
    //     PrimaryVariables values;

    //     values[Indices::pressureIdx] = 1e5;
    //     values[Indices::saturationIdx] = 0.0;

    //     if (onUpperBoundary(globalPos))
    //     {
    //         if (dropletSolver_->isCoupledWithDroplet(globalPos))
    //         {

    //             const auto& droplet = dropletSolver_->droplet();
    //             values[Indices::pressureIdx] = 1e5 + dropletSolver_->Pc(droplet);
    //             values[Indices::saturationIdx] = 1.0;
    //         }
    //         else
    //         {
    //             values[Indices::pressureIdx] = 1e5;
    //             values[Indices::saturationIdx] = 0.0;
    //         }
    //     }

    //     return values;
    // }

    PrimaryVariables dirichlet(const Element &element, const SubControlVolume &scv) const
    {
        PrimaryVariables values;

        const auto& globalPos = scv.dofPosition();

        values[Indices::pressureIdx] = 1e5;
        values[Indices::saturationIdx] = 0.0;

        if (onUpperBoundary(globalPos))
        {
            if (dropletSolver_->isCoupledWithDroplet(globalPos))
            {

                const auto& droplet = dropletSolver_->droplet();
                values[Indices::pressureIdx] = 1e5 + dropletSolver_->Pc(droplet);
                values[Indices::saturationIdx] = 1.0;
            }
            else
            {
                values[Indices::pressureIdx] = 1e5;
                values[Indices::saturationIdx] = 0.0;
            }
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    // NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    // {
    //     NumEqVector values(0.0);
    //     if (onUpperBoundary(globalPos))
    //     {
    //         if (dropletSolver_->isCoupledWithDroplet(globalPos))
    //         {
    //             values[Indices::conti0EqIdx] = -0.002; // kg/(m*s)
    //         }
    //     }
    //     return values;
    // }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        // forward it to the interface with only the global position
        // const auto& globalPos = scvf.ipGlobal();
        NumEqVector values(0.0);
        // if (onUpperBoundary(globalPos))
        // {
        //     if (dropletSolver_->isCoupledWithDroplet(globalPos))
        //     {
        //         // Scalar dropletInfiltrationRate = dropletSolver_->dropletMassFlux(element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache); //kg/s
        //         // dropletInfiltrationRate /= scvf.area();
        //         values[Indices::conti0EqIdx] = dropletInfiltrationRate; // kg/(m*s) or kg/(m2*s)
        //         // std::cout<<"------dropletInfiltrationRate-------"<<dropletInfiltrationRate<<std::endl;
        //     }
        // }
        return values;
    }

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = 1e5;
        values[Indices::saturationIdx] = 0.0;


        // if (onUpperBoundary(globalPos))
        // {
        //     if (dropletSolver_->isCoupledWithDroplet(globalPos))
        //     {

        //         const auto& droplet = dropletSolver_->droplet();
        //         values[Indices::pressureIdx] = 1e5 + dropletSolver_->Pc(droplet);
        //         values[Indices::saturationIdx] = 1.0;
        //     }
        // }

        return values;
    }

    void setDropSolver(std::shared_ptr<DropSolver> dropletSolver)
    { dropletSolver_ = dropletSolver; }

    auto dropSolver()
    { return dropletSolver_; }

    auto dropSolver() const
    { return dropletSolver_; }


    bool onInlet(const GlobalPosition &globalPos) const
    {
        Scalar width = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        Scalar lambda = (this->gridGeometry().bBoxMax()[0] - globalPos[0])/width;
        return onUpperBoundary(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }


   bool onLeftBoundary(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
    }

    bool onCircleBoundary(const GlobalPosition &globalPos) const
    {
        return std::hypot(globalPos[0]-tabletCenter_[0], globalPos[2]-tabletCenter_[2]) > tabletRadius_ - eps_;
    }

    GlobalPosition tabletCenter() const
    { return tabletCenter_; }

    GlobalPosition tabletCenterUpperBoundary() const
    {
        auto center = tabletCenter_;
        center[1] = this->gridGeometry().bBoxMax()[1];
        return center;
    }

    Scalar tabletRadius() const
    { return tabletRadius_; }
private:
    static constexpr Scalar eps_ = 1e-6;
    const GlobalPosition tabletCenter_;
    const Scalar tabletRadius_;


    std::shared_ptr<DropSolver> dropletSolver_;
};

} // end namespace Dumux

#endif
