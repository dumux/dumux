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

namespace Dumux {
template <class TypeTag>
class DarcySubProblem;

namespace Properties {

// create new property for projection operator
template<class TypeTag, class MyTypeTag>
struct MortarProjector { using type = UndefinedProperty; };

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
    using Projector = GetPropType<TypeTag, Properties::MortarProjector>;

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
    {
        const auto mortarVariable = getParamFromGroup<std::string>("Mortar", "VariableType");
        if (mortarVariable == "Pressure")
            useDirichletAtInterface_ = true;
        else if (mortarVariable == "Flux")
            useDirichletAtInterface_ = false;

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
        BoundaryTypes values;

        values.setAllNeumann();
        if (!isOnNegativeMortarSide_ && onLowerBoundary_(scv.dofPosition()))
            values.setAllDirichlet();
        if (isOnNegativeMortarSide_ && onUpperBoundary_(scv.dofPosition()))
            values.setAllDirichlet();

        if (useDirichletAtInterface_)
            if (isOnMortarInterface(scv.dofPosition()))
                values.setAllDirichlet();

        return values;
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
        BoundaryTypes values;

        values.setAllNeumann();
        if (!isOnNegativeMortarSide_ && onLowerBoundary_(scvf.ipGlobal()))
            values.setAllDirichlet();
        if (isOnNegativeMortarSide_ && onUpperBoundary_(scvf.ipGlobal()))
            values.setAllDirichlet();

        if (useDirichletAtInterface_)
            if (isOnMortarInterface(scvf.ipGlobal()))
                values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Evaluate the source term at a given position.
     */
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        static const auto topSource = getParam<Scalar>("Problem.TopSource");
        const auto xMax = this->fvGridGeometry().bBoxMax()[0];
        if (isOnNegativeMortarSide_ && !useHomogeneousSetup_)
            return NumEqVector(topSource*(xMax - globalPos[0]));
        return NumEqVector(0.0);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet sub-control volume face
     * \note This overload is for the box scheme
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolume& scv) const
    {
        if (useDirichletAtInterface_ && isOnMortarInterface(scv.dofPosition()))
            return PrimaryVariables(computeMortarPressure_(element));
        return dirichletAtPos(scv.dofPosition());
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet sub-control volume face
     * \note This overload is for cell-centered schemes
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        if (useDirichletAtInterface_ && isOnMortarInterface(scvf.ipGlobal()))
            return PrimaryVariables(computeMortarPressure_(element));
        return dirichletAtPos(scvf.ipGlobal());
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet segment
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        static const auto topBoundaryPressure = getParam<Scalar>("Problem.TopBoundaryPressure");
        static const auto bottomBoundaryPressure = getParam<Scalar>("Problem.BottomBoundaryPressure");
        if (!useHomogeneousSetup_)
        {
            if (!isOnNegativeMortarSide_ && onLowerBoundary_(globalPos))
                return PrimaryVariables(bottomBoundaryPressure);
            else if (isOnNegativeMortarSide_ && onUpperBoundary_(globalPos))
                return  PrimaryVariables(topBoundaryPressure);
        }

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
        if ( isOnMortarInterface(scvf.ipGlobal()) )
        {
            auto flux = mortarProjector_().integrateMortarVariable(element);
            // turn it into flux per area
            flux /= mortarProjector_().overlapArea(element);
            // scale with density (mortar variable is velocity)
            flux *= elemVolVars[scvf.insideScvIdx()].density();

            if (isOnNegativeMortarSide_)
                flux *= -1.0;

            return NumEqVector(flux);
        }
        else
            return NumEqVector(0.0);
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


    //! set the pointer to the projector class
    void setMortarProjector(std::shared_ptr<const Projector> p)
    {
        projector_ = p;
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

private:
    //! computes the mortar pressure (if pressure coupling is used)
    Scalar computeMortarPressure_(const Element& element) const
    {
        assert(mortarProjector_().hasOverlapWithMortar(element));
        return mortarProjector_().integrateMortarVariable(element)
              / mortarProjector_().overlapArea(element);
    }

    //! Returns true if position is on lower domain boundary
    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] < this->fvGridGeometry().bBoxMin()[dimWorld-1] + eps_; }

    //! Returns true if position is on upper domain boundary
    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps_; }

    //! Returns a const reference to the projector class
    const Projector& mortarProjector_() const
    {
        if (!projector_)
            DUNE_THROW(Dune::InvalidStateException, "Projector pointer not set");
        return *projector_;
    }

    Scalar eps_;
    std::string problemName_;
    std::shared_ptr<const Projector> projector_;

    bool isOnNegativeMortarSide_;
    bool useHomogeneousSetup_;
    bool useDirichletAtInterface_;
};
} // end namespace Dumux

#endif //DUMUX_DARCY_SUBPROBLEM_HH
