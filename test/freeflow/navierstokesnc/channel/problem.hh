// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup NavierStokesNCTests
 * \brief Channel flow test for the multi-component staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_CHANNEL_NC_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_NC_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/massproblem.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesNCTests
 * \brief  Test problem for the one-phase (Navier-)Stokes model.
 *
 * Flow from left to right in a channel is considered. A parabolic velocity
 * profile is set at the left boundary, while the pressure is set to
 * a fixed value on the right boundary. The top and bottom boundaries
 * represent solid walls with no-slip/no-flow conditions.
 */
template <class TypeTag>
class ChannelNCTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
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
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;


    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    static constexpr bool useMoles = ModelTraits::useMoles();
    static constexpr bool isNonisothermal = ModelTraits::enableEnergyBalance();
    static constexpr auto compIdx = 1;

public:
    ChannelNCTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , eps_(1e-6)
    {
        FluidSystem::init();
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
    }

    VelocityVector faceVelocity(const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const SubControlVolumeFace& scvf) const
    {
        const Scalar velocity = oldStaggeredVelocity_[scvf.index()];
        auto faceNormal = scvf.unitOuterNormal();
        for (int i=0; i<faceNormal.size(); ++i)
            faceNormal[i] = std::abs(faceNormal[i]);
        return velocity * faceNormal;
    }

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is substracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return 1.1e5; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if (!isInlet_(globalPos))
            values.setAllNeumann();

        if (isInlet_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
            values.setDirichlet(Indices::conti0EqIdx + compIdx);
            if constexpr (isNonisothermal)
                values.setDirichlet(Indices::temperatureIdx);
        }

        if (isOutlet_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
            values.setDirichlet(Indices::conti0EqIdx + compIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    DirichletValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto& globalPos = scvf.ipGlobal();
        DirichletValues values = initialAtPos(globalPos);

        // give the system some time so that the pressure can equilibrate, then start the injection of the tracer
        if (isInlet_(globalPos))
        {
            //values[Indices::pressureIdx] = this->couplingManager().cellPressure(element, scvf);
            values[Indices::pressureIdx] = oldStaggeredMass_[scvf.insideScvIdx()][0];

            if (time() >= 10.0 || inletVelocity_  < eps_)
            {
                Scalar moleFracTransportedComp = 1e-3;
                if constexpr (useMoles)
                    values[Indices::conti0EqIdx+compIdx] = moleFracTransportedComp;
                else
                {
                    Scalar averageMolarMassPhase = moleFracTransportedComp * FluidSystem::molarMass(compIdx)
                                               + (1. - moleFracTransportedComp)  * FluidSystem::molarMass(1-compIdx);
                    values[Indices::conti0EqIdx+compIdx] = moleFracTransportedComp * FluidSystem::molarMass(compIdx)
                                            / averageMolarMassPhase;
                }

            if constexpr (isNonisothermal)
                values[Indices::temperatureIdx] = 293.15;
            }
        }
        //if (isOutlet_(globalPos))
        //    values[Indices::pressureIdx] = 1.1e+5;

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);

        const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();
        values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf) * insideDensity * scvf.unitOuterNormal();

        using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;
        if (isOutlet_(scvf.ipGlobal()))
            values = FluxHelper::scalarOutflowFlux( *this, element, fvGeometry, scvf, elemVolVars);

        return values;
    }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    InitialValues initialAtPos(const GlobalPosition& globalPos) const
    {
        InitialValues values(0.0);

        values[Indices::pressureIdx] = 1.1e+5;
        if constexpr (isNonisothermal)
            values[Indices::temperatureIdx] = 283.15;

        return values;
    }

    /*!
     * \brief A parabolic velocity profile.
     *
     * \param y The position where the velocity is evaluated.
     * \param vMax The profile's maxmium velocity.
     */
    Scalar parabolicProfile(const Scalar y, const Scalar vMax) const
    {
        const Scalar yMin = this->gridGeometry().bBoxMin()[1];
        const Scalar yMax = this->gridGeometry().bBoxMax()[1];
        return  vMax * (y - yMin)*(yMax - y) / (0.25*(yMax - yMin)*(yMax - yMin));
    }

    void setTimeLoop(TimeLoopPtr timeLoop)
    { timeLoop_ = timeLoop; }

    Scalar time() const
    { return timeLoop_->time(); }

    //! Enable internal Dirichlet constraints
    static constexpr bool enableInternalDirichletConstraints()
    { return false; }

    /*!
     * \brief Tag a degree of freedom to carry internal Dirichlet constraints.
     *        If true is returned for a dof, the equation for this dof is replaced
     *        by the constraint that its primary variable values must match the
     *        user-defined values obtained from the function internalDirichlet(),
     *        which must be defined in the problem.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     */
    std::bitset<DirichletValues::dimension> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<DirichletValues::dimension> values;
        for (const auto& intersection : intersections(this->gridGeometry().gridView(), element))
        {
            const auto center = intersection.geometry().center();
            if (intersection.boundary() && center[0]
                    < this->gridGeometry().bBoxMin()[0] + eps_)
                values.set(0);
        }
        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    DirichletValues internalDirichlet(const Element& element, const SubControlVolume& scv) const
    //{ return DirichletValues(1.1e+5); }
    { return DirichletValues(oldStaggeredMass_[scv.dofIndex()][0]); }

    void setOldStaggered(Dune::FieldVector<Scalar, 5000> oldStaggeredVelocity,
            Dune::FieldVector<Dune::FieldVector<Scalar, 2>, 1250> oldStaggeredMass)
    {
        oldStaggeredVelocity_ = oldStaggeredVelocity;
        oldStaggeredMass_ = oldStaggeredMass;
    }

private:
    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    const Scalar eps_;
    Scalar inletVelocity_;
    TimeLoopPtr timeLoop_;
    Dune::FieldVector<Scalar, 5000> oldStaggeredVelocity_;
    Dune::FieldVector<Dune::FieldVector<Scalar, 2>, 1250> oldStaggeredMass_;
};
} // end namespace Dumux

#endif
