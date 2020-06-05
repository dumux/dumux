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
 * \ingroup GeomechanicsTests
 * \brief Definition of a test problem for the poro-elastic model.
 */

#ifndef DUMUX_POROELASTIC_PROBLEM_HH
#define DUMUX_POROELASTIC_PROBLEM_HH

#include <dune/common/fmatrix.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/discretization/box.hh>
#include <dumux/geomechanics/poroelastic/model.hh>
#include <dumux/geomechanics/fvproblem.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class PoroElasticProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct TestPoroElastic { using InheritsFrom = std::tuple<PoroElastic, BoxModel>; };
} // end namespace TTag
// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TestPoroElastic> { using type = Dune::YaspGrid<2>; };
// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TestPoroElastic> { using type = Dumux::PoroElasticProblem<TypeTag>; };
// The fluid phase consists of one constant component
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TestPoroElastic>
{
    using type = Dumux::FluidSystems::OnePLiquid< GetPropType<TypeTag, Properties::Scalar>,
                                                  Dumux::Components::Constant<0, GetPropType<TypeTag, Properties::Scalar>> >;
};
// The spatial parameters property
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TestPoroElastic>
{
    using type = PoroElasticSpatialParams< GetPropType<TypeTag, Properties::Scalar>,
                                           GetPropType<TypeTag, Properties::GridGeometry> >;
};
} // end namespace Properties

/*!
 * \ingroup GeomechanicsTests
 * \brief Problem definition for the deformation of a poro-elastic body.
 */
template<class TypeTag>
class PoroElasticProblem : public GeomechanicsFVProblem<TypeTag>
{
    using ParentType = GeomechanicsFVProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr Scalar pi = M_PI;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GradU = Dune::FieldMatrix<Scalar, dim, dimWorld>;

public:
    PoroElasticProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    //! Returns the temperature in the domain.
    static constexpr Scalar temperature()
    { return 273.15; }

    //! Evaluates the initial value for a control volume.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! Evaluates the boundary conditions for a Dirichlet boundary segment.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Returns the effective fluid density.
     *
     * \param globalPos The global position
     */
    Scalar effectiveFluidDensityAtPos(const GlobalPosition& globalPos) const
    {
        // This test uses the constant component, obtain density only once
        using FS = GetPropType<TypeTag, Properties::FluidSystem>;
        static const Scalar rho = FS::density( effectivePorePressureAtPos(globalPos), temperature() );
        return rho;
    }

    /*!
     * \brief Returns the effective pore pressure
     *
     * \note We use the x-displacement as pressure solution. The shift to
     *       higher values is done to see a mor pronounced effect in stresses.
     *
     * \param globalPos The global position
     */
    Scalar effectivePorePressureAtPos(const GlobalPosition& globalPos) const
    { return exactSolution(globalPos)[0] + 10; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume.
     */
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        using std::sin;
        using std::cos;

        const auto ipGlobal = scv.center();
        const auto x = ipGlobal[0];
        const auto y = ipGlobal[1];

        // the lame parameters (we know they only depend on position here)
        const auto& lameParams = this->spatialParams().lameParamsAtPos(scv.center());
        const auto lambda = lameParams.lambda();
        const auto mu = lameParams.mu();

        // precalculated products
        const Scalar pi_2 = 2.0*pi;
        const Scalar pi_2_square = pi_2*pi_2;
        const Scalar cos_2pix = cos(pi_2*x);
        const Scalar sin_2pix = sin(pi_2*x);
        const Scalar cos_2piy = cos(pi_2*y);
        const Scalar sin_2piy = sin(pi_2*y);

        const Scalar dE11_dx = -2.0*sin_2piy;
        const Scalar dE22_dx = pi_2_square*cos_2pix*cos_2piy;
        const Scalar dE11_dy = pi_2*(1.0-2.0*x)*cos_2piy;
        const Scalar dE22_dy = -1.0*pi_2_square*sin_2pix*sin_2piy;
        const Scalar dE12_dy = 0.5*pi_2_square*(cos_2pix*cos_2piy - (x-x*x)*sin_2piy);
        const Scalar dE21_dx = 0.5*((1.0-2*x)*pi_2*cos_2piy - pi_2_square*sin_2pix*sin_2piy);

        // compute exact divergence of sigma
        PrimaryVariables divSigma(0.0);
        divSigma[Indices::momentum(/*x-dir*/0)] = lambda*(dE11_dx + dE22_dx) + 2*mu*(dE11_dx + dE12_dy);
        divSigma[Indices::momentum(/*y-dir*/1)] = lambda*(dE11_dy + dE22_dy) + 2*mu*(dE21_dx + dE22_dy);
        return divSigma;
    }

    /*!
     * \brief Evaluates the exact displacement to this problem at a given position.
     */
    PrimaryVariables exactSolution(const GlobalPosition& globalPos) const
    {
        using std::sin;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        PrimaryVariables exact(0.0);
        exact[Indices::momentum(/*x-dir*/0)] = (x-x*x)*sin(2*pi*y);
        exact[Indices::momentum(/*y-dir*/1)] = sin(2*pi*x)*sin(2*pi*y);
        return exact;
    }

    /*!
     * \brief Evaluates the exact displacement gradient to this problem at a given position.
     */
    GradU exactGradient(const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::cos;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        static constexpr int xIdx = Indices::momentum(/*x-dir*/0);
        static constexpr int yIdx = Indices::momentum(/*y-dir*/1);

        GradU exactGrad(0.0);
        exactGrad[xIdx][xIdx] = (1-2*x)*sin(2*pi*y);
        exactGrad[xIdx][yIdx] = (x - x*x)*2*pi*cos(2*pi*y);
        exactGrad[yIdx][xIdx] = 2*pi*cos(2*pi*x)*sin(2*pi*y);
        exactGrad[yIdx][yIdx] = 2*pi*sin(2*pi*x)*cos(2*pi*y);
        return exactGrad;
    }

private:
    static constexpr Scalar eps_ = 3e-6;
};

} // end namespace Dumux

#endif
