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
 * \ingroup NavierStokesTests
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model.
 *
 * The channel is either modeled in 3D or in 2D, using an additional wall friction term
 * to mimic the 3D behavior of the flow.
 */

#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROBLEM_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/freeflow/navierstokes/problem.hh>

namespace Dumux {

/*!
 * \brief Test problem for the one-phase (Navier-) Stokes model in a 3D or pseudo 3D channel.
 *
 * Flow from left to right in a three-dimensional channel is considered. At the inlet (left)
 * and outlet (right) fixed values for pressure are set.
 * The channel is confined by solid walls at all other sides of the domain which corresponds
 * to no-slip/no-flow conditions.
 * The value of an analytical solution for the given flow configuration is furthermore provided.
 * For sake of efficiency, the 3D problem can be reduced to a two-dimensional one by including
 * an additional wall friction term to the momentum balance (Flekkoy et al., 1995 \cite flekkoy1995a).
 */
template <class TypeTag>
class ThreeDChannelTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = typename ParentType::NumEqVector;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    template<class GridData>
    ThreeDChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager, GridData gridData)
    : ParentType(gridGeometry, couplingManager)
    {
        deltaP_ = getParam<Scalar>("Problem.DeltaP");

        boundaryFlag_.resize(gridGeometry->numScvf(), 0);
        auto fvGeometry = localView(*gridGeometry);
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& scvf : scvfs(fvGeometry))
                if (scvf.boundary())
                    boundaryFlag_[scvf.index()]
                        = gridData->getBoundaryDomainMarker(scvf.boundaryFlag());
        }
    }

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
    { return 273.15 + 10; } // 10C

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

        if constexpr (ParentType::isMomentumProblem())
        {
            values.setAllDirichlet();
            if (isOutlet_(scvf) || isInlet_(scvf))
                values.setAllNeumann();
        }
        else
            values.setNeumann(Indices::conti0EqIdx);

        return values;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        static_assert(!ParentType::isMomentumProblem());
        values.setNeumann(Indices::conti0EqIdx);
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // no-flow/no-slip
        return PrimaryVariables(0.0);
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
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            if (isOutlet_(scvf) || isInlet_(scvf))
            {
                const auto& fluxVarCache = elemFluxVarsCache[scvf];
                using Tensor = Dune::FieldMatrix<Scalar, dim>;
                Tensor gradV(0.0);
                for (const auto& scv : scvs(fvGeometry))
                {
                    const auto& volVars = elemVolVars[scv];
                    for (int dir = 0; dir < dim; ++dir)
                        gradV[dir].axpy(volVars.velocity(dir), fluxVarCache.gradN(scv.indexInElement()));
                }

                static const bool enableUnsymmetrizedVelocityGradient
                    = getParamFromGroup<bool>(this->paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

                values = enableUnsymmetrizedVelocityGradient ?
                    NumEqVector(0.0)//mv(gradV, scvf.unitOuterNormal())
                    : mv(getTransposed(gradV), scvf.unitOuterNormal()); //mv(gradV + getTransposed(gradV), scvf.unitOuterNormal());

                values *= -this->effectiveViscosity(element, fvGeometry, scvf);

                const auto referencePressure = this->referencePressure(element, fvGeometry, scvf);
                const auto p = (isInlet_(scvf) ? deltaP_ : 0.0) - referencePressure;
                values.axpy(p, scvf.unitOuterNormal());
            }
        }
        else
        {
            if (isInlet_(scvf) || isOutlet_(scvf))
            {
                const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();
                values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf) * insideDensity * scvf.unitOuterNormal();
            }
        }

        return values;
    }

    // \}

private:

    bool isInlet_(const SubControlVolumeFace& scvf) const
    { return boundaryFlag_[scvf.index()] == 1; }

    bool isOutlet_(const SubControlVolumeFace& scvf) const
    { return boundaryFlag_[scvf.index()] == 2; }

    static constexpr Scalar eps_ = 1e-10;
    Scalar deltaP_;

    std::vector<int> boundaryFlag_;
};

} // end namespace Dumux

#endif
