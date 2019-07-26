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
 * \brief A simple Darcy test problem (cell-centered finite volume method).
 */
#ifndef DUMUX_DARCY_SUBPROBLEM_HH
#define DUMUX_DARCY_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "spatialparams_darcy.hh"

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "mortarvariabletype.hh"

namespace Dumux {
template <class TypeTag>
class DarcySubProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct DarcyOneP { using InheritsFrom = std::tuple<OneP>; };
struct DarcyOnePTpfa { using InheritsFrom = std::tuple<DarcyOneP, CCTpfaModel>; };
struct DarcyOnePMpfa { using InheritsFrom = std::tuple<DarcyOneP, CCMpfaModel>; };
struct DarcyOnePBox { using InheritsFrom = std::tuple<DarcyOneP, BoxModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOneP> { using type = Dumux::DarcySubProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyOneP>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<0, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyOneP>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<Scalar, 2>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyOneP>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePDarcySpatialParams<FVGridGeometry, Scalar>;
};
} // end namespace Properties

template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr bool isBox = FVGridGeometry::discMethod == DiscretizationMethod::box;

public:
    //! export spatial parameters type
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    DarcySubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                    std::shared_ptr<SpatialParams> spatialParams,
                    const std::string& paramGroup)
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , eps_(1e-7)
    , isOnNegativeMortarSide_(getParamFromGroup<bool>(paramGroup, "Problem.IsOnNegativeMortarSide"))
    , useHomogeneousSetup_(false)
    {
        const auto mv = getParamFromGroup<std::string>("Mortar", "VariableType");
        mortarVariableType_ = mv == "Pressure" ? OnePMortarVariableType::pressure : OnePMortarVariableType::flux;
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
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
      * \brief Specifies which kind of boundary condition should be
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scvf The boundary sub control volume face
      */
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolume &scv) const
    {
        if (isOnMortarInterface(scv.dofPosition()))
        {
            BoundaryTypes values;
            if (mortarVariableType_ == OnePMortarVariableType::pressure)
                values.setAllDirichlet();
            else
                values.setAllNeumann();
            return values;
        }
        else
            return boundaryTypesAtPos(scv.dofPosition());
    }

    /*!
      * \brief Specifies which kind of boundary condition should be
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scvf The boundary sub control volume face
      */
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolumeFace &scvf) const
    {
        if (isOnMortarInterface(scvf.ipGlobal()))
        {
            BoundaryTypes values;
            if (mortarVariableType_ == OnePMortarVariableType::pressure)
                values.setAllDirichlet();
            else
                values.setAllNeumann();
            return values;
        }
        else
            return boundaryTypesAtPos(scvf.ipGlobal());
    }

    /*!
      * \brief Specifies which kind of boundary condition should be
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scvf The boundary sub control volume face
      */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluate the exact solution at a given position.
     */
    NumEqVector exact(const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::cos;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        return NumEqVector( sin(M_PI*x)*cos(M_PI*(y-0.5)) + x*x );
    }

    /*!
     * \brief Evaluate the exact normal velocity at a given position.
     */
    GlobalPosition exactFlux(const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::cos;

        const auto x = globalPos[0];
        const auto y = globalPos[1];
        const auto& k = this->spatialParams().permeabilityAtPos(globalPos);

        return GlobalPosition( {-1.0*k*M_PI*cos(M_PI*x)*cos(M_PI*(y-0.5)) - 2.0*x*k, k*M_PI*sin(M_PI*x)*sin(M_PI*(y-0.5))} );
    }

    /*!
     * \brief Evaluate the source term at a given position.
     */
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        if (!useHomogeneousSetup_)
        {
            using std::sin;
            using std::cos;

            const auto x = globalPos[0];
            const auto y = globalPos[1];
            const auto& k = this->spatialParams().permeabilityAtPos(globalPos);

            return NumEqVector( k*M_PI*M_PI*(sin(M_PI*x)*cos(M_PI*(y-0.5)) - 2.0*k + sin(M_PI*x)*cos(M_PI*(y-0.5)))  );
        }

        return NumEqVector(0.0);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet sub-control volume face
     * \note This overload is for the box scheme
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolume& scv) const
    {
        assert(mortarProjection_.size() == this->fvGridGeometry().numDofs());
        if (isOnMortarInterface(scv.dofPosition()))
            return PrimaryVariables(mortarProjection_[scv.dofIndex()]);
        return dirichletAtPos(scv.dofPosition());
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet sub-control volume face
     * \note This overload is for cell-centered schemes
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        assert(mortarProjection_.size() == this->fvGridGeometry().numDofs());
        if (isOnMortarInterface(scvf.ipGlobal()))
            return PrimaryVariables( mortarProjection_[scvf.insideScvIdx()] );
        return dirichletAtPos(scvf.ipGlobal());
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet segment
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        if (!useHomogeneousSetup_)
            return exact(globalPos);
        else
            return PrimaryVariables(0.0);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scvf The boundary sub control volume face
     *
     * For this method, the \a values variable stores primary variables.
     */
    template<class ElementVolumeVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
    {
        const auto globalPos = scvf.ipGlobal();
        if ( isOnMortarInterface(globalPos) )
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());

            auto flux = mortarProjection_[insideScv.elementIndex()];
            flux *= elemVolVars[insideScv].density();

            return isOnNegativeMortarSide_ ? NumEqVector(-1.0*flux) : NumEqVector(flux);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "This test should have no Neumann BCS");
    }

    // \}

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param element The element
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }


    //! set the projected mortar solution
    void setMortarProjection(SolutionVector p)
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
        return (isOnNegativeMortarSide_ && onBottomBoundary_(globalPos))
               || (!isOnNegativeMortarSide_ && onTopBoundary_(globalPos));
    }

    //! Returns true if this domain is on the "negative" side of mortar
    bool isOnNegativeMortarSide() const
    { return isOnNegativeMortarSide_; }

    //! Define the meaning of the mortar variable
    void setMortarVariableType(OnePMortarVariableType mv)
    {
        mortarVariableType_ = mv;
    }

private:
    //! Returns true if position is on lower domain boundary
    bool onBottomBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] < this->fvGridGeometry().bBoxMin()[dimWorld-1] + eps_; }

    //! Returns true if position is on upper domain boundary
    bool onTopBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps_; }

    //! Returns true if position is on right domain boundary
    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_; }

    //! Returns true if position is on left domain boundary
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_; }

    Scalar eps_;
    std::string problemName_;
    SolutionVector mortarProjection_;

    bool isOnNegativeMortarSide_;
    bool useHomogeneousSetup_;
    OnePMortarVariableType mortarVariableType_;
};

} // end namespace Dumux

#endif //DUMUX_DARCY_SUBPROBLEM_HH
