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
 * \ingroup PoromechanicsTests
 * \brief The poro-elastic sub-problem in the el2p coupled problem.
 */
#ifndef DUMUX_POROELASTIC_SUBPROBLEM_HH
#define DUMUX_POROELASTIC_SUBPROBLEM_HH

#include <dune/common/fmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>

#include <dumux/geomechanics/fvproblem.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The poro-elastic sub-problem in the el2p coupled problem.
 */
template<class TypeTag>
class PoroElasticSubProblem : public GeomechanicsFVProblem<TypeTag>
{
    using ParentType = GeomechanicsFVProblem<TypeTag>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GradU = Dune::FieldMatrix<Scalar, dim, dimWorld>;

    using StressType = GetPropType<TypeTag, Properties::StressType>;
    using StressTensor = typename StressType::StressTensor;
    using ForceVector = typename StressType::ForceVector;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using FluxVarCache = typename
    GridVariables::GridFluxVariablesCache::FluxVariablesCache;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
public:
    PoroElasticSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                          std::shared_ptr<GetPropType<TypeTag, Properties::SpatialParams>> spatialParams,
                          const std::string& paramGroup = "PoroElastic")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        stress_[0].resize(gridGeometry->gridView().size(0));
        stress_[1].resize(gridGeometry->gridView().size(0));
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    //! Evaluates the initial value for a control volume.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! Evaluates the boundary conditions for a Dirichlet boundary segment.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! Evaluates the boundary conditions for a Neumann boundary segment.
    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables force(0.0); // [Pa]
        force[1] = - 2260 * 9.8 * -globalPos[1];
        force[0] = 0.6 * force[1];
        //std::cout << "At pos"<< globalPos << " Neumann value: " << force << std::endl;
        return force;
    }
    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if (globalPos[0] < eps_)
            values.setDirichlet(0);

        if (globalPos[1] < -2500 + eps_)
            values.setDirichlet(1);
        return values;
    }

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume.
     */
    PrimaryVariables source(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    { return PrimaryVariables(0.0); }

    void calculateStress(const GridVariables& gv,
                         const SolutionVector& sol)
    {
        auto& gg = this->gridGeometry();
        for (const auto& element : elements(gg.gridView()))
        {
            auto fvGeometry = localView(gg);
            auto elemVolVars = localView(gv.curGridVolVars());

            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, sol);

            // evaluate flux variables cache at cell center
            FluxVarCache fluxVarCache;
            fluxVarCache.update(*this, element, fvGeometry, elemVolVars, element.geometry().center());

            const auto sigma = StressType::stressTensor(*this, element, fvGeometry, elemVolVars, fluxVarCache);

            const auto eIdx = gg.elementMapper().index(element);
            for(int dir = 0; dir < dim; dir++)
            {
                stress_[dir][eIdx] = sigma[dir];
            }
        }
    }

    auto& getStress(const int& dir) const
    {
        return stress_[dir];
    }

private:
    static constexpr Scalar eps_ = 3e-6;
    std::string problemName_;
    std::array<std::vector<ForceVector>,2> stress_;
};

} // end namespace Dumux

#endif
