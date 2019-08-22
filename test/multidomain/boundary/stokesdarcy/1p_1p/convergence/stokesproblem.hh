// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup NavierStokesTests
 * \brief A simple Stokes test problem for the staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>

namespace Dumux
{
template <class TypeTag>
class StokesSubProblem;

namespace Properties
{
NEW_TYPE_TAG(StokesOneP, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokes));

// the fluid system
SET_PROP(StokesOneP, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
SET_TYPE_PROP(StokesOneP, Grid, Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(StokesOneP, Problem, Dumux::StokesSubProblem<TypeTag> );

SET_BOOL_PROP(StokesOneP, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(StokesOneP, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(StokesOneP, EnableGridVolumeVariablesCache, true);
}

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem.
 *
 * Horizontal flow from left to right with a parabolic velocity profile.
 */
template <class TypeTag>
class StokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

public:
    StokesSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "Stokes"), eps_(1e-6), couplingManager_(couplingManager)
    {
        createAnalyticalSolution_();
    }

   /*!
     * \name Problem parameters
     */
    // \{

    bool shouldWriteRestartFile() const
    { return false; }

   /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C

   /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector source(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::cos;
        using std::sin;

        source[Indices::conti0EqIdx]  = -sin(2*M_PI*x);
        source[Indices::momentumXBalanceIdx] = -2*M_PI*y*y*sin(4*M_PI*x)
                                               -2*y*sin(2*M_PI*x)
                                               +2*M_PI*cos(2*M_PI*x) ;
        source[Indices::momentumYBalanceIdx] =  -2*M_PI*y*y*cos(2*M_PI*x)
                                                -4*M_PI*M_PI*sin(2*M_PI*x);

        return source;
    }
    // \}

   /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;

        const auto& globalPos = scvf.dofPosition();

        values.setDirichlet(Indices::velocityXIdx);
        values.setDirichlet(Indices::velocityYIdx);

        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values.setCouplingNeumann(Indices::conti0EqIdx);
            values.setCouplingNeumann(Indices::momentumYBalanceIdx);
            values.setBJS(Indices::momentumXBalanceIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element
     * \param scvf The sub control volume face
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        return analyticalSolution(globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeomentry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values[Indices::conti0EqIdx] = couplingManager().couplingData().massCouplingCondition(fvGeometry, elemVolVars, elemFaceVars, scvf);
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(fvGeometry, elemVolVars, elemFaceVars, scvf);
        }

        return values;
    }

    // \}

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return analyticalSolution(globalPos);
    }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar permeability(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().couplingData().darcyPermeability(scvf);
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().problem(CouplingManager::darcyIdx).spatialParams().beaversJosephCoeffAtPos(scvf.center());
    }

    // \}

    /*!
      * \brief Returns the analytical solution for the pressure
      */
     auto& getAnalyticalPressureSolution() const
     {
         return analyticalPressure_;
     }

    /*!
      * \brief Returns the analytical solution for the velocity
      */
     auto& getAnalyticalVelocitySolution() const
     {
         return analyticalVelocity_;
     }

     /*!
       * \brief Returns the analytical solution for the velocity at the faces
       */
      auto& getAnalyticalVelocitySolutionOnFace() const
      {
          return analyticalVelocityOnFace_;
      }

private:
    /*!
    * \brief Return the analytical solution of the problem at a given position
    *
    * \param globalPos The global position
    */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = exactPressure(globalPos);

        auto velocity = exactVelocity(globalPos);
        values[Indices::velocityXIdx] = velocity[0];
        values[Indices::velocityYIdx] = velocity[1];

        return values;
    }

    Scalar exactPressure(const GlobalPosition &globalPos) const
    {
        using std::sin;
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return -y*y*sin(2*M_PI*x)*sin(2*M_PI*x);
    }

    GlobalPosition exactVelocity(const GlobalPosition &globalPos) const
    {
        using std::sin;
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        GlobalPosition velocity(0.0);
        velocity[0] = y;
        velocity[1] = -y*sin(2*M_PI*x);
        return velocity;
    }

    /*!
      * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
      */
     void createAnalyticalSolution_()
     {
         analyticalPressure_.resize(this->fvGridGeometry().numCellCenterDofs());
         analyticalVelocity_.resize(this->fvGridGeometry().numCellCenterDofs());
         analyticalVelocityOnFace_.resize(this->fvGridGeometry().numFaceDofs());

         for (const auto& element : elements(this->fvGridGeometry().gridView()))
         {
             auto fvGeometry = localView(this->fvGridGeometry());
             fvGeometry.bindElement(element);
             for (auto&& scv : scvs(fvGeometry))
             {
                 auto ccDofIdx = scv.dofIndex();
                 auto ccDofPosition = scv.dofPosition();
                 auto analyticalSolutionAtCc = analyticalSolution(ccDofPosition);

                 // velocities on faces
                 for (auto&& scvf : scvfs(fvGeometry))
                 {
                     const auto faceDofIdx = scvf.dofIndex();
                     const auto faceDofPosition = scvf.center();
                     const auto dirIdx = scvf.directionIndex();
                     const auto analyticalSolutionAtFace = analyticalSolution(faceDofPosition);
                     analyticalVelocityOnFace_[faceDofIdx][dirIdx] = analyticalSolutionAtFace[Indices::velocity(dirIdx)];
                 }

                 analyticalPressure_[ccDofIdx] = analyticalSolutionAtCc[Indices::pressureIdx];

                 for(int dirIdx = 0; dirIdx < ModelTraits::dim(); ++dirIdx)
                     analyticalVelocity_[ccDofIdx][dirIdx] = analyticalSolutionAtCc[Indices::velocity(dirIdx)];
             }
         }
      }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_; }

    Scalar eps_;

    std::shared_ptr<CouplingManager> couplingManager_;
    std::vector<Scalar> analyticalPressure_;
    std::vector<GlobalPosition> analyticalVelocity_;
    std::vector<GlobalPosition> analyticalVelocityOnFace_;
};
} //end namespace

#endif // DUMUX_STOKES_SUBPROBLEM_HH
