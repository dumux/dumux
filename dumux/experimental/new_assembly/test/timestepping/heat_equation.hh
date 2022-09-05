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
 * \brief Types required for solving the heat equation using TPFA.
 */
#ifndef DUMUX_TEST_TIMESTEPPING_HEAT_EQUATION_HH
#define DUMUX_TEST_TIMESTEPPING_HEAT_EQUATION_HH

#include <memory>
#include <cmath>

#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>
#include <dumux/experimental/new_assembly/dumux/assembly/operatorweights.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/gridvariables.hh>

namespace Dumux {

template<typename Scalar>
struct HeatEquationProblem
{
    template<typename Position>
    bool isDirichlet(const Position& p) const
    { return p[0] < 1e-6 || p[0] > 1.0 - 1e-6; }

    template<typename Position>
    Scalar dirichletValue(const Position& p) const
    { return p[0]; }

    Scalar heatConductivity() const
    { return 1.0; }

    Scalar heatCapacity() const
    { return 1.0; }
};

// Custom grid variables to define a local view
template<typename GridGeometry,
         typename DofVector,
         typename IndexStrategy>
class HeatEquationGridVariables
: public GridVariables<GridGeometry, DofVector, IndexStrategy>
{
    using ParentType = GridVariables<GridGeometry, DofVector, IndexStrategy>;

public:
    using Scalar = Dumux::LinearSystem::ScalarType<DofVector>;

    class LocalView
    {
    public:
        explicit LocalView(const HeatEquationGridVariables& gridVars)
        : gridVars_(gridVars)
        {}

        template<typename GGLocalView>
        void bind(const GGLocalView) &
        {}

        template<typename GGLocalView>
        LocalView bind(const GGLocalView) &&
        { return std::move(*this); }

        const HeatEquationGridVariables& gridVariables() const
        { return gridVars_; }

    private:
        const HeatEquationGridVariables& gridVars_;
    };

    template<typename Initializer>
    HeatEquationGridVariables(std::shared_ptr<const GridGeometry> gridGeometry,
                              std::shared_ptr<const IndexStrategy> indexStrategy,
                              const Initializer& initializer)
    : ParentType(gridGeometry, indexStrategy, initializer)
    {}

    const HeatEquationProblem<Scalar>& problem() const
    { return problem_; }

private:
    HeatEquationProblem<Scalar> problem_;
};

// custom local assembler for our heat equation
template<typename GG, typename GV>
class HeatEquationLocalAssembler
{
    using LocalGridGeometry = typename GG::LocalView;
    using LocalGridVariables = typename GV::LocalView;
    using Scalar = typename GV::Scalar;

    using SubControlVolume = typename GG::SubControlVolume;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;

public:
    using GridGeometry = GG;
    using GridVariables = GV;
    using OperatorWeights = Dumux::OperatorWeights<Scalar>;

    HeatEquationLocalAssembler(const LocalGridGeometry& localGeom,
                               const LocalGridVariables& localVars,
                               const OperatorWeights& weights)
    : localGeom_(localGeom)
    , localVars_(localVars)
    , weights_(weights)
    {}

    template<typename Matrix>
    void addJacobianEntries(Matrix& jacobian) const
    { addJacobianEntries_(jacobian); }

    template<typename Vector>
    void addResidualEntries(Vector& residual) const
    { addResidualEntries_(residual); }

    template<typename Matrix, typename Vector>
    void addJacobianAndResidualEntries(Matrix& jacobian, Vector& residual) const
    {
        addJacobianEntries(jacobian);
        addResidualEntries(residual);
    }

private:
    decltype(auto) problem_() const
    { return localVars_.gridVariables().problem(); }

    bool isNeumannFace_(const SubControlVolumeFace& scvf) const
    { return scvf.boundary() && !problem_().isDirichlet(scvf.ipGlobal()); }

    template<typename MultiIndex>
    decltype(auto) getDofIndex_(const MultiIndex& i) const
    { return localVars_.gridVariables().getDofIndex(i); }

    template<typename Vector>
    void addResidualEntries_(Vector& residual) const
    {
        for (const auto& scv : scvs(localGeom_))
            LinearSystem::add(
                residual,
                getDofIndex_(MultiIndex{scv.dofIndex(), 0}),
                evalTemporal_(scv)
            );

        for (const auto& scvf : scvfs(localGeom_))
            LinearSystem::add(
                residual,
                getDofIndex_(MultiIndex{scvf.insideScvIdx(), 0}),
                evalSpatial_(scvf)
            );
    }

    Scalar evalTemporal_(const SubControlVolume& scv) const
    { return storageDerivative_(scv)*getDofValue_(scv); }

    Scalar evalSpatial_(const SubControlVolumeFace& scvf) const
    {
        if (weights_.spatialWeight)
            return weights_.spatialWeight.value()*computeFaceFlux_(scvf);
        return 0.0;
    }

    Scalar storageDerivative_(const SubControlVolume& scv) const
    {
        const auto cp = problem_().heatCapacity();
        if (weights_.temporalWeight)
            return weights_.temporalWeight.value()*cp*scv.volume();
        return 0.0;
    }

    template<typename Matrix>
    void addJacobianEntries_(Matrix& jacobian) const
    {
        for (const auto& scv : scvs(localGeom_))
            LinearSystem::add(
                jacobian,
                getDofIndex_(MultiIndex{scv.dofIndex(), 0}),
                getDofIndex_(MultiIndex{scv.dofIndex(), 0}),
                storageDerivative_(scv)
            );

        for (const auto& scvf : scvfs(localGeom_))
        {
            if (isNeumannFace_(scvf))
                continue;

            const auto weight = weights_.spatialWeight
                                ? weights_.spatialWeight.value()
                                : 0.0;
            const auto deriv = weight*computeFaceTransmissibility_(scvf);
            LinearSystem::add(
                jacobian,
                getDofIndex_(MultiIndex{scvf.insideScvIdx(), 0}),
                getDofIndex_(MultiIndex{scvf.insideScvIdx(), 0}),
                deriv
            );

            if (!scvf.boundary())
                LinearSystem::add(
                    jacobian,
                    getDofIndex_(MultiIndex{scvf.insideScvIdx(), 0}),
                    getDofIndex_(MultiIndex{scvf.outsideScvIdx(), 0}),
                    -1.0*deriv
                );
        }
    }

    double computeFaceFlux_(const SubControlVolumeFace& scvf) const
    {
        if (isNeumannFace_(scvf))
            return 0.0;

        const auto tij = computeFaceTransmissibility_(scvf);
        const auto ui = getInsideDofValue_(scvf);
        const auto uj = scvf.boundary() ? problem_().dirichletValue(scvf.ipGlobal())
                                        : getOutsideDofValue_(scvf);
        return tij*(ui - uj);
    }

    double computeFaceTransmissibility_(const SubControlVolumeFace& scvf) const
    {
        const auto heatConductivity = problem_().heatConductivity();
        const auto& insideScv = localGeom_.scv(scvf.insideScvIdx());
        const auto insideT = Dumux::computeTpfaTransmissibility(
            scvf, insideScv, heatConductivity, 1.0
        )*scvf.area();
        if (scvf.boundary())
            return insideT;

        const auto& outsideScv = localGeom_.scv(scvf.outsideScvIdx());
        const auto outsideT = Dumux::computeTpfaTransmissibility(
            scvf, outsideScv, heatConductivity, 1.0
        )*scvf.area();
        return insideT*outsideT/(outsideT - insideT);
    }

    double getInsideDofValue_(const SubControlVolumeFace& scvf) const
    { return getDofValue_(localGeom_.scv(scvf.insideScvIdx())); }

    double getOutsideDofValue_(const SubControlVolumeFace& scvf) const
    { return getDofValue_(localGeom_.scv(scvf.outsideScvIdx())); }

    double getDofValue_(const SubControlVolume& scv) const
    {
        return LinearSystem::get(
            localVars_.gridVariables().dofs(),
            getDofIndex_(MultiIndex{scv.dofIndex(), 0})
        );
    }

    const LocalGridGeometry& localGeom_;
    const LocalGridVariables& localVars_;
    OperatorWeights weights_;
};

} // namespace Dumux

#endif
