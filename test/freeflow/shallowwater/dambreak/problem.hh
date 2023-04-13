// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief A test for the Shallow water model (wet dam break).
 */
#ifndef DUMUX_DAM_BREAK_TEST_PROBLEM_HH
#define DUMUX_DAM_BREAK_TEST_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/shallowwater/problem.hh>
#include <dumux/flux/shallowwater/riemannproblem.hh>
#include <dumux/flux/shallowwater/exactriemann.hh>

namespace Dumux {

/*!
 * \ingroup Shallow water equations model
 * \ingroup ImplicitTestProblems
 *
 * \brief A simple dam break test (1D wet dam break).
 *
 * The domain is a long rectangle with a gate in the middle. On the left
 * side the water depth is larger than on the right side.
 * All boundaries are set to no-flow.
 *
 * This problem uses the \ref ShallowWaterModel.
 *
 */
template <class TypeTag>
class DamBreakProblem : public ShallowWaterProblem<TypeTag>
{
    using ParentType = ShallowWaterProblem<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using NeumannFluxes = Dumux::NumEqVector<PrimaryVariables>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    DamBreakProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
        exactWaterDepth_.resize(gridGeometry->numDofs(), 0.0);
        exactVelocityX_.resize(gridGeometry->numDofs(), 0.0);
    }

    //! Get the analytical water depth
    const std::vector<Scalar>& getExactWaterDepth()
    {
        return exactWaterDepth_;
    }

    //! Get the analytical velocity
    const std::vector<Scalar>& getExactVelocityX()
    {
        return exactVelocityX_;
    }

    //! Update the analytical solution
    template<class SolutionVector, class GridVariables>
    void updateAnalyticalSolution(const SolutionVector& curSol,
                                  const GridVariables& gridVariables,
                                  const Scalar time)
    {
        //compute solution for all elements
        auto fvGeometry = localView(this->gridGeometry());
        auto elemVolVars = localView(gridVariables.curGridVolVars());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, curSol);

            const auto& globalPos = element.geometry().center();

            //compute the position s and the initial water depth at the gate, velocities are zero
            const Scalar s = (globalPos[0] - gatePosition_)/time;
            const Scalar waterDepthLeft =  initialWaterDepthLeft_;
            const Scalar waterDepthRight =  initialWaterDepthRight_;
            const auto gravity = this->spatialParams().gravity(globalPos);

            auto riemannResult = ShallowWater::exactRiemann(waterDepthLeft,
                                                            waterDepthRight,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            gravity,
                                                            s);

            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            exactWaterDepth_[eIdx] = riemannResult.waterDepth;
            exactVelocityX_[eIdx] = riemannResult.velocityX;
        }
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }


    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the boundary type is set
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Specifies the neumann boundary
     * \param element
     * \param fvGeometry
     * \param elemVolVars
     * \param elemFluxVarsCache
     * \param scvf
     */
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const ElementFluxVariablesCache& elemFluxVarsCache,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);

        // we need the Riemann invariants to compute the values depending of the boundary type
        // since we use a weak imposition we do not have a dirichlet value. We impose fluxes
        // based on q,h, etc. computed with the Riemann invariants
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& nxy = scvf.unitOuterNormal();
        const auto gravity = this->spatialParams().gravity(scvf.center());

        auto riemannFlux = ShallowWater::riemannProblem(insideVolVars.waterDepth(),
                                                        insideVolVars.waterDepth(),
                                                        insideVolVars.velocity(0),
                                                        -insideVolVars.velocity(0),
                                                        insideVolVars.velocity(1),
                                                        -insideVolVars.velocity(1),
                                                        insideVolVars.bedSurface(),
                                                        insideVolVars.bedSurface(),
                                                        gravity,
                                                        nxy);

        values[Indices::massBalanceIdx] = riemannFlux[0];
        values[Indices::velocityXIdx]   = riemannFlux[1];
        values[Indices::velocityYIdx]   = riemannFlux[2];

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial values for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * \param globalPos The position for which the boundary type is set
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {

        PrimaryVariables values(0.0);

        values[0] = initialWaterDepthRight_;
        values[1] = 0.0;
        values[2] = 0.0;

        // water level on the left side of the gate
        if (globalPos[0] < 10.0 + eps_)
        {
            values[0] = initialWaterDepthLeft_;
        }

        //water level on the right side of the gate
        else
        {
            values[0] = initialWaterDepthRight_;
        }

        return values;
    };

    // \}

private:

    std::vector<Scalar> exactWaterDepth_;
    std::vector<Scalar> exactVelocityX_;

    static constexpr Scalar initialWaterDepthLeft_ = 4.0;
    static constexpr Scalar initialWaterDepthRight_ = 1.0;
    static constexpr Scalar channelLenght_ = 20.0;
    static constexpr Scalar gatePosition_ = 10.0;

    static constexpr Scalar eps_ = 1.0e-6;
    std::string name_;
};

} //end namespace Dumux

#endif
