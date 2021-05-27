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
 * \ingroup BoundaryTests
 * \brief A simple Stokes test problem for the staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

#include <utility>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/geometry/makegeometry.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>

#include <test/freeflow/navierstokes/l2error.hh>

#include "mortarvariabletype.hh"

namespace Dumux {
template <class TypeTag>
class StokesSubProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct StokesOneP { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StokesOneP>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<0, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesOneP> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesOneP> { using type = Dumux::StokesSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup BoundaryTests
 * \brief Test problem for the one-phase (Navier-) Stokes problem.
 *
 * Horizontal flow from left to right with a parabolic velocity profile.
 */
template <class TypeTag>
class StokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    // extract sub-vector for cell-centered values (used for mortar data)
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using CellSolutionVector = std::decay_t<decltype( std::declval<SolutionVector>()[FVGridGeometry::cellCenterIdx()] )>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<ModelTraits::numEq()>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    StokesSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, const std::string& paramGroup)
    : ParentType(fvGridGeometry, paramGroup)
    , eps_(1e-6)
    , isOnNegativeMortarSide_(getParamFromGroup<bool>(paramGroup, "Problem.IsOnNegativeMortarSide"))
    {
        const auto mv = getParamFromGroup<std::string>("Mortar", "VariableType");
        mortarVariableType_ = mv == "Pressure" ? OnePMortarVariableType::pressure
                                               : OnePMortarVariableType::flux;

        problemName_  =  getParamFromGroup<std::string>(paramGroup, "Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    { return problemName_; }

   /*!
     * \name Problem parameters
     */
    // \{

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C

    //! compute the l2 error of the solution
    template <class SolutionVector>
    PrimaryVariables calculateL2Error(const SolutionVector& curSol) const
    {
        using L2Error = NavierStokesTestL2Error<Scalar, ModelTraits, PrimaryVariables>;
        const auto l2error = L2Error::calculateL2Error(*this, curSol);
        const int numCellCenterDofs = this->gridGeometry().numCellCenterDofs();
        const int numFaceDofs = this->gridGeometry().numFaceDofs();
        std::cout << std::setprecision(8) << "** L2 error (abs) for "
                << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
                << std::scientific
                << "L2 errors = " << l2error.first[Indices::pressureIdx] << "   "
                << l2error.first[Indices::velocityXIdx] << "   "
                << l2error.first[Indices::velocityYIdx]
                << std::endl;

        return l2error.first;
    }

    // \}

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume face.
     */
    using ParentType::source;
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const ElementFaceVariables& elemFaceVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);

        if (!useHomogeneousSetup_)
        {
            using std::sin;
            const auto& quad = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(element.geometry().type(), 5);
            for (auto&& qp : quad)
            {
                const auto pos = element.geometry().global(qp.position());
                Scalar x = pos[0];
                Scalar y = pos[1];
                auto integrationElement = element.geometry().integrationElement(qp.position());
                // source[Indices::conti0EqIdx]  += -1.0*sin(2.0*M_PI*x) * qp.weight()*integrationElement;
                source[Indices::conti0EqIdx]  += 0.0*qp.weight()*integrationElement;
            }
            source /= scv.volume();
        }

        return source;
    }

    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector source(const Element &element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace &scvf) const
    {
        NumEqVector source(0.0);

        if (!useHomogeneousSetup_)
        {
            using std::cos;
            using std::sin;

            std::vector<Dune::FieldVector<Scalar, 3>> corners(4);
            corners[0] = {scvf.corner(0)[0], scvf.corner(0)[1], 0};
            corners[1] = {scvf.corner(1)[0], scvf.corner(1)[1], 0};
            corners[2] = corners[0];
            corners[3] = corners[1];
            corners[2][scvf.directionIndex()] = element.geometry().center()[scvf.directionIndex()];
            corners[3][scvf.directionIndex()] = element.geometry().center()[scvf.directionIndex()];

            const auto geometry = makeDuneQuadrilaterial(corners);
            const auto& quad = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(geometry.type(), 5);

            for (auto&& qp : quad)
            {
                Dune::FieldVector<Scalar, 3> globalPos = geometry.global(qp.position());
                Scalar x = globalPos[0];
                Scalar y = globalPos[1];

                auto integrationElement = geometry.integrationElement(qp.position());
                // source[Indices::momentumXBalanceIdx] += (-4.0*M_PI*y*y*sin(2.0*M_PI*x)*cos(2.0*M_PI*x)
                //                                         -2.0*y*sin(2.0*M_PI*x) + 2.0*M_PI*cos(2.0*M_PI*x))
                //                                          * qp.weight()*integrationElement;
                //
                // source[Indices::momentumYBalanceIdx] += (-2.0*M_PI*y*y*cos(2.0*M_PI*x) - 2.0*M_PI*2.0*M_PI*y*sin(2.0*M_PI*x))
                //                                           * qp.weight()*integrationElement;
                source[Indices::momentumXBalanceIdx] += 0.0*qp.weight()*integrationElement;
                source[Indices::momentumYBalanceIdx] += 0.0*qp.weight()*integrationElement;
             }

            source /= 0.5*element.geometry().volume();
        }

        return source;
    }

    Scalar exactPressure(const GlobalPosition &globalPos) const
    {
        using std::sin;
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        // return -y*y*sin(2.0*M_PI*x)*sin(2.0*M_PI*x);
        return 2.0*x + y - 1;
    }

    GlobalPosition exactVelocity(const GlobalPosition &globalPos) const
    {
        using std::sin;
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        GlobalPosition velocity(0.0);
        // velocity[0] = y;
        // velocity[1] = -y*sin(2.0*M_PI*x);
        velocity[0] = (y-1.0)*(y-1.0) + x*(y-1.0) + 3.0*x - 1;
        velocity[1] = x*(x-1.0) - 0.5*(y-1.0)*(y-1.0) - 3.0*y + 1;
        return velocity;
    }

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

        if (isOnMortarInterface(globalPos))
        {
            if (mortarVariableType_ == OnePMortarVariableType::flux)
                values.setDirichlet(scvf.directionIndex());
            else
                values.setNeumann(scvf.directionIndex());
            values.setBeaversJoseph(1 - scvf.directionIndex());
        }
        else
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluates the velocity boundary conditions for a Dirichlet sub-control volume face.
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        if (mortarVariableType_ == OnePMortarVariableType::flux && isOnMortarInterface(scvf.ipGlobal()))
        {
            PrimaryVariables values(0.0);
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            values[scvf.directionIndex()] = mortarProjection_[eIdx];
            return values;
        }

        if (!useHomogeneousSetup_)
            return analyticalSolution(scvf.ipGlobal());
        else
            return PrimaryVariables(0.0);
    }

    /*!
     * \brief Define pressure level inside a cell
     */
    bool isDirichletCell(const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const SubControlVolume& scv,
                         int pvIdx) const
    {
        if (mortarVariableType_ == OnePMortarVariableType::flux
            && pvIdx == 2 && scv.dofIndex() == this->gridGeometry().gridView().size(0)-1)
                return true;
        return false;

    }

    /*!
     * \brief Evaluates the pressure boundary conditions for a sub-control volume
     *        touchgin a dirichlet boundary segment.
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolume& scv) const
    { return analyticalSolution(scv.center()); }

    /*!
     * \brief Define Dirichlet boundary conditions on a position on a boundary segment.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos); }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
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

        if (isOnMortarInterface(scvf.ipGlobal()) && mortarVariableType_ != OnePMortarVariableType::flux)
        {
            // apply mortar pressure to momentum balance
            assert(mortarProjection_.size() == this->gridGeometry().gridView().size(0));
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            values[scvf.directionIndex()] = mortarProjection_[eIdx]*scvf.directionSign();

            // mass flux
            const auto v = elemFaceVars[scvf].velocitySelf()*scvf.directionSign();
            values[Indices::conti0EqIdx] = v*elemVolVars[scvf.insideScvIdx()].density();
        }

        return values;
    }

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluates the initial value for a control volume.
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! returns the intrinsic permeability
    Scalar permeability(const Element& element, const SubControlVolumeFace& scvf) const
    { return 1.0; }

    //! returns the alpha value
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        static const Scalar alpha = getParam<Scalar>("SpatialParams.AlphaBeaversJoseph");
        return alpha;
    }

    //! set the pointer to the projector class
    void setMortarProjection(CellSolutionVector p)
    {
        mortarProjection_ = p;
    }

    //! Set whether or not the homogeneous system is solved
    void setUseHomogeneousSetup(bool value)
    {
        useHomogeneousSetup_ = value;
    }

    //! Returns true if a position if on the mortar interface
    bool isOnMortarInterface(const GlobalPosition& globalPos) const
    {
        return (isOnNegativeMortarSide_ && onLowerBoundary_(globalPos))
               || (!isOnNegativeMortarSide_ && onUpperBoundary_(globalPos));
    }

    //! Returns true if this domain is on the "negative" side of mortar
    bool isOnNegativeMortarSide() const
    { return isOnNegativeMortarSide_; }

    //! Define the meaning of the mortar variable
    void setMortarVariableType(OnePMortarVariableType mv)
    {
        mortarVariableType_ = mv;
    }

    CellSolutionVector mortarProjection() const
    { return mortarProjection_; }

    /*!
    * \brief Return the analytical solution of the problem at a given position
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

private:
    bool onLeftBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    Scalar eps_;
    std::string problemName_;
    CellSolutionVector mortarProjection_;

    bool isOnNegativeMortarSide_;
    bool useHomogeneousSetup_;
    OnePMortarVariableType mortarVariableType_;
};
} // end namespace Dumux

#endif // DUMUX_STOKES_SUBPROBLEM_HH
