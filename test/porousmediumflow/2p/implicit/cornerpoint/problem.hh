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
 * \ingroup TwoPTests
 * \brief The properties for the 2p cornerpoint test
 */
#ifndef DUMUX_TWOP_CORNERPOINT_TEST_PROBLEM_HH
#define DUMUX_TWOP_CORNERPOINT_TEST_PROBLEM_HH

#include <opm/grid/CpGrid.hpp>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/2p/incompressiblelocalresidual.hh>

#include "spatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class TwoPCornerPointTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct TwoPCornerPoint { using InheritsFrom = std::tuple<CCTpfaModel, TwoP>; };
} // end namespace TTag

// Set the grid type
SET_TYPE_PROP(TwoPCornerPoint, Grid, Dune::CpGrid);

// Set the problem type
SET_TYPE_PROP(TwoPCornerPoint, Problem, TwoPCornerPointTestProblem<TypeTag>);

// the local residual containing the analytic derivative methods
SET_TYPE_PROP(TwoPCornerPoint, LocalResidual, TwoPIncompressibleLocalResidual<TypeTag>);

// Set the fluid system
SET_PROP(TwoPCornerPoint, FluidSystem)
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Set the spatial parameters
SET_PROP(TwoPCornerPoint, SpatialParams)
{
private:
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = TwoPCornerPointTestSpatialParams<FVGridGeometry, Scalar>;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TwoPCornerPoint> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TwoPCornerPoint> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::TwoPCornerPoint> { static constexpr bool value = false; };
} // end namespace Properties

/*!
 * \ingroup TwoPTests
 * \brief The incompressible 2p cornerpoint test problem.
 */
template<class TypeTag>
class TwoPCornerPointTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum { dimWorld = GridView::dimensionworld };

public:
    TwoPCornerPointTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                               std::shared_ptr<typename ParentType::SpatialParams> spatialParams)
    : ParentType(fvGridGeometry, spatialParams)
    {
        gravity_ = {0, 0, 9.81};
        injectionElement_ = getParam<int>("Problem.InjectionElement");
        injectionRate_ = getParam<Scalar>("Problem.InjectionRate");
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes bcTypes;

        // set no-flux on top and bottom, hydrostatic on the rest
        // use the scvf normal to decide
        const auto& normal = scvf.unitOuterNormal();
        using std::abs;
        if (abs(normal[dimWorld-1]) > 0.5)
            bcTypes.setAllNeumann();
        else
            bcTypes.setAllDirichlet();

        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initialAtPos(globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    //! \copydoc Dumux::FVProblem::source()
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        NumEqVector values(0.0);

        int eIdx = this->fvGridGeometry().gridView().indexSet().index(element);
        if (eIdx == injectionElement_)
            values[FluidSystem::phase1Idx] = injectionRate_/element.geometry().volume();

        return values;
    }

    const GlobalPosition& gravity() const
    {
        return gravity_;
    }

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        // hydrostatic pressure
        Scalar densityW = 1000;
        values[Indices::pressureIdx] = 1e5 + densityW*(this->gravity()*globalPos);
        values[Indices::saturationIdx] = 0.0;

        return values;
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    {
        return 293.15;
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    template<class VTKWriter>
    void addFieldsToWriter(VTKWriter& vtk)
    {
        const auto numElements = this->fvGridGeometry().gridView().size(0);

        permX_.resize(numElements);
        permZ_.resize(numElements);

        vtk.addField(permX_, "PERMX [mD]");
        vtk.addField(permZ_, "PERMZ [mD]");

        const auto& gridView = this->fvGridGeometry().gridView();
        for (const auto& element : elements(gridView))
        {
            const auto eIdx = this->fvGridGeometry().elementMapper().index(element);

            // transfer output to mD = 9.86923e-16 m^2
            permX_[eIdx] = this->spatialParams().permeabilityX(eIdx)/9.86923e-16;
            permZ_[eIdx] = this->spatialParams().permeabilityZ(eIdx)/9.86923e-16;
        }
    }

private:
    GlobalPosition gravity_;
    int injectionElement_;
    Scalar injectionRate_;
    std::vector<Scalar> permX_, permZ_;
};

} // end namespace Dumux

#endif
