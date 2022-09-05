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
 * \ingroup Assembly
 * \copydoc Dumux::CCFVLocalAssembler
 */
#ifndef DUMUX_ASSEMBLY_CCFV_LOCAL_ASSEMBLER_HH
#define DUMUX_ASSEMBLY_CCFV_LOCAL_ASSEMBLER_HH

#include <optional>
#include <concepts>

#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/variables.hh>
#include <dumux/experimental/new_assembly/dumux/assembly/operatorweights.hh>

namespace Dumux {
namespace Concepts {

template<typename T, typename Scalar = double>
concept NumEqVector = Indexable<T> and ScaleAssignable<Scalar>;

} // namespace Concepts

/*!
 * \ingroup Assembly
 * \brief Local assembler for PDEs discretized with cell-centered schemes.
 */
template<typename GridGeometry, typename GridVariables>
class CCFVLocalAssembler
{
    using LocalGridGeometry = typename GridGeometry::LocalView;
    using LocalGridVariables = typename GridVariables::LocalView;

    using Scalar = Variables::ScalarType<GridVariables>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    static constexpr int numEq = GridVariables::Model::numEq;
    using NumEqVector  = Dune::FieldVector<Scalar, numEq>;
    using FluxesVector = Dune::ReservedVector<
        NumEqVector,
        LocalGridGeometry::maxElementStencilSize
    >;

public:
    using OperatorWeights = Dumux::OperatorWeights<Scalar>;
    using ElementResidualVector = Dune::ReservedVector<
        NumEqVector,
        LocalGridGeometry::maxNumElementScvs
    >;

    CCFVLocalAssembler(const LocalGridGeometry& localGeom,
                       LocalGridVariables& localVars,
                       std::optional<OperatorWeights> weights = {})
    : localGeom_(localGeom)
    , localVars_(localVars)
    , weights_(weights)
    {}

    template<LinearSystem::Interoperable Residual>
    void addResidualEntries(Residual& residual) const
    { addScvResidualEntries_(residual, evalLocalResidual()); }

    template<LinearSystem::Interoperable Matrix>
    void addJacobianEntries(Matrix& jacobian)
    { addJacobianEntries_(jacobian, evalLocalResidual()); }

    template<LinearSystem::Interoperable Matrix,
             LinearSystem::Interoperable Vector>
    void addJacobianAndResidualEntries(Matrix& jacobian, Vector& residual)
    {
        const auto origLocalResidual = evalLocalResidual();
        addJacobianEntries_(jacobian, origLocalResidual)
        addScvResidualEntries_(residual, origLocalResidual);
    }

    ElementResidualVector evalLocalResidual() const
    {
        ElementResidualVector result(localGeom_.numScv());
        evalLocalResidual_(result);
        return result;
    }

private:
    const auto& gridVariables_() const
    { return localVars_.gridVariables(); }

    const auto& model_() const
    { return gridVariables_().model(); }

    const auto& problem_() const
    { return gridVariables_().problem(); }

    template<LinearSystem::Interoperable Residual>
    void addScvResidualEntries_(Residual& res,
                                const ElementResidualVector& values) const
    {
        for (const auto& scv : scvs(localGeom_))
            addResidualEntries_(res, values[scv.localDofIndex()], scv.dofIndex());
    }

    template<LinearSystem::Interoperable Residual>
    void addResidualEntries_(Residual& res,
                             const NumEqVector& values,
                             const std::size_t dofIdx) const
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            LinearSystem::add(
                res,
                gridVariables_().getDofIndex(MultiIndex{dofIdx, eqIdx}),
                values[eqIdx]
            );
    }

    template<LinearSystem::Interoperable Matrix>
    void addJacobianEntries_(Matrix& jacobian,
                             const ElementResidualVector& origResidual)
    {
        FluxesVector origNeighborFluxes;
        if (doSpatial_())
            origNeighborFluxes = computeNeighborScvfFluxes_();
        auto neighborFluxDerivs = origNeighborFluxes;

        auto derivs = origResidual;
        auto elemVarsView = deflectableView(localVars_);

        for (const auto& scvJ : scvs(localGeom_))
        {
            const auto globalJ = scvJ.dofIndex();
            for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
            {
                const auto eps = 1e-8;
                const auto origPrivar = localVars_.getPrimaryVariable(scv, pvIdx);
                const auto deflectedPriVar = origPrivar + eps;

                elemVarsView.deflect(scvJ, pvIdx, deflectedPriVar);

                evalLocalResidual_(derivs);
                derivs -= origResidual;
                derivs /= eps;

                for (const auto& scvI : scvs(localGeom_))
                {
                    const auto globalI = scvI.dofIndex();
                    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                        LinearSystem::add(
                            jacobian,
                            gridVariables_().getDofIndex(MultiIndex{globalI, eqIdx}),
                            gridVariables_().getDofIndex(MultiIndex{globalJ, pvIdx}),
                            derivs[scvI.localDofIndex()][0]
                        );
                }

                if (doSpatial_())
                {
                    evalNeighborScvfFluxes_(neighborFluxDerivs);
                    neighborFluxDerivs -= origNeighborFluxes;
                    neighborFluxDerivs /= eps;

                    int localNeighborScvfIdx = 0;
                    for (const auto& scvfI : neighborScvfs(localGeom_))
                    {
                        const auto scvI = localGeom_.insideScv(scvfI);
                        const auto globalI = scvI.dofIndex();
                        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                            LinearSystem::add(
                                jacobian,
                                gridVariables_().getDofIndex(MultiIndex{globalI, eqIdx}),
                                gridVariables_().getDofIndex(MultiIndex{globalJ, pvIdx}),
                                neighborFluxDerivs[scvI.localDofIndex()][localNeighborScvfIdx++]
                            );
                    }
                }

                // TODO: Restore always? Nex update may be enough?
                // Also, how to tell that only vol vars should be deflected,
                // e.g. when spatial = false?
                elemVarsView.restore();
            }
        }
    }

    auto evalSources_(const SubControlVolume& scv) const
    {
        auto sources = problem_.source(localGeom_, localVars_, scv);
        sources *= scv.volume()*spatialWeight_();
        return sources;
    }

    auto evalStorage_(const SubControlVolume& scv) const
    {
        auto storage = evalOperator_(model_().storageOperator(), scv);
        storage *= scv.volume()*temporalWeight_();
        return storage;
    }

    auto evalFluxes_(const SubControlVolumeFace& scvf) const
    {
        const auto scaleFactor = scvf.area()*spatialWeight_();
        if (isNeumannFace_(scvf))
        {
            auto values = problem_().neumann(localGeom_, localVars_, scvf);
            values *= scaleFactor;
            return values;
        }
        else
        {
            auto values = evalOperator_(model_().fluxOperator(), scvf);
            values *= scaleFactor;
            return values;
        }
    }

    FluxesVector computeNeighborScvfFluxes_() const
    {
        FluxesVector result;
        evalNeighborScvfFluxes_(result);
        return result;
    }

    void evalNeighborScvfFluxes_(FluxesVector& neighborFluxes) const
    {
        int i = 0;
        neighborFluxes.resize(localGeom_.numNeighborScvfs());
        for (const auto& scvfJ : neighborScvfs(localGeom_))
            neighborFluxes[i++] = evalFluxes_(scvfJ);
    }

    void evalLocalResidual_(ElementResidualVector& res) const
    {
        res = 0.0;
        if (doSpatial_())
        {
            for (const auto& scv : scvs(localGeom_))
                result[scv.localDofIndex()] += evalSources_(scv);

            for (const auto& scvf : scvfs(localGeom_))
                result[scv.localDofIndex()] += evalFluxes_(scvf);
        }

        if (doTemporal_())
            for (const auto& scv : scvs(localGeom_))
                result[scv.localDofIndex()] += evalStorage_(scv);
    }

    template<typename Operator, typename SCE> requires(
        std::invocable<Operator, const LG&, const LV&, const SCE&> and
        std::same_as<NumEqVector, std::invoke_result_t<Operator, const LG&, const LV&, const SCE&>>)
    NumEqVector evalOperator_(const Operator& op,
                              const SCE& subControlEntity) const
    { return op(localGeom_, localVars_, subControlEntity); }

    bool isNeumannFace_(const SubControlVolume& scvf) const
    { return problem_().boundaryTypes(localGeom_, scvf).hasNeumann(); }

    bool doSpatial_() const
    { return weights_ && weights_.spatialWeight; }

    bool doTemporal_() const
    { return weights_ && weights_.temporalWeight; }

    Scalar temporalWeight_() const
    { return doTemporal_() ? weights_.temporalWeight : 1.0; }

    Scalar spatialWeight_() const
    { return doSpatial_() ? weights_.spatialWeight : 1.0; }

    const ElementGeometry& localGeom_;
    ElementVariables& localVars_;
    std::optional<OperatorWeights> weights_;
};

} // namespace Dumux

#endif
