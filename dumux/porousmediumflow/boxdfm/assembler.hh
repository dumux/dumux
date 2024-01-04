// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDFMModel
 * \ingroup Assembly
 * \copydoc Dumux::BoxDfmAssembler
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_BOX_DFM_ASSEMBLER_HH
#define DUMUX_POROUS_MEDIUM_FLOW_BOX_DFM_ASSEMBLER_HH

#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/common/exceptions.hh>

#include <dumux/geometry/volume.hh>
#include <dumux/common/numericdifferentiation.hh>

#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/volvardeflectionhelper_.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup BoxDFMModel
 * \ingroup Assembly
 * \brief A linear system assembler for the box-dfm model, which allows to incorporate
 *        constraints across barriers.
 * \note This assembler is only required if barriers are present. The standard Dumux::FVAssembler
 *       works fine for simulations without barriers.
 * \note For now, this assembler is restricted to numeric differentation and implicit assembly.
 */
template<class TypeTag, class BarrierFlux>
class BoxDfmAssembler : public FVAssembler<TypeTag, DiffMethod::numeric, /*implicit*/true>
{
    using ParentType = FVAssembler<TypeTag, DiffMethod::numeric, true>;
    static constexpr int dim = GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension;

public:
    using typename ParentType::GridVariables;
    using typename ParentType::GridGeometry;
    using typename ParentType::JacobianMatrix;
    using typename ParentType::SolutionVector;
    using typename ParentType::ResidualType;
    using typename ParentType::Scalar;

private:
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    struct BarrierIntegrationPoint
    {
        const GlobalPosition& integrationPoint;
        const GlobalPosition& normalVector;
        const typename GlobalPosition::value_type weight;

        const FVElementGeometry& insideFVGeometry;
        const FVElementGeometry& outsideFVGeometry;

        const SubControlVolume& insideScv;
        const SubControlVolume& outsideScv;

        ElementVolumeVariables& insideElemVolVars;
        ElementVolumeVariables& outsideElemVolVars;
    };

    template<class... Args>
    BoxDfmAssembler(BarrierFlux&& barrierFlux, Args&&... args)
    : ParentType(std::forward<Args>(args)...)
    , barrierFlux_(std::move(barrierFlux))
    {}

    //! Assembles the global Jacobian of the residual and the residual for the current solution.
    template<class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(const SolutionVector& curSol, const PartialReassembler* partialReassembler = nullptr)
    {
        ParentType::assembleJacobianAndResidual(curSol, partialReassembler);
        assembleBarrierEntries_(curSol, this->jacobian(), this->residual());
    }

    //! Assembles only the global Jacobian of the residual.
    void assembleJacobian(const SolutionVector& curSol)
    {
        ParentType::assembleJacobian(curSol);
        DUNE_THROW(Dune::NotImplemented, "assembleJacobian(x)");
    }

    //! compute the residuals using the internal residual
    void assembleResidual(const SolutionVector& curSol)
    {
        ParentType::assembleResidual(curSol);
        DUNE_THROW(Dune::NotImplemented, "assembleResidual(x)");
    }

    //! assemble a residual r
    void assembleResidual(ResidualType& r, const SolutionVector& curSol) const
    {
        ParentType::assembleResidual(r, curSol);
        DUNE_THROW(Dune::NotImplemented, "assembleResidual(r, x)");
    }

private:
    BarrierFlux barrierFlux_;

    void assembleBarrierEntries_(const SolutionVector& curSol, JacobianMatrix& A, ResidualType& r) const
    {
        static constexpr auto numEq = SolutionVector::value_type::size();

        const auto computeBarrierFlux = [&] (const BarrierIntegrationPoint& ip) {
            auto flux = barrierFlux_(ip);
            flux *= ip.weight;
            return flux;
        };

        visitBarriers_(curSol, [&] (BarrierIntegrationPoint ip) {
            const bool isOnBoundary = this->gridGeometry().dofOnBoundary(ip.insideScv.dofIndex());
            const auto bcTypes = this->problem().boundaryTypes(ip.insideFVGeometry.element(), ip.insideScv);

            const auto flux = computeBarrierFlux(ip);
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                if (!bcTypes.isDirichlet(eqIdx) || !isOnBoundary)
                    r[ip.insideScv.dofIndex()][eqIdx] += flux[eqIdx];

            const auto addDerivatives = [&] (auto& elemSol, auto& evv, const auto& fvGeo) {
                auto partialDerivs = flux;
                auto deflectionHelper = Detail::VolVarsDeflectionHelper{
                    [&] (const auto& scv) -> VolumeVariables& { return evv[scv]; },
                    fvGeo,
                    true
                };

                for (const auto& scv : scvs(fvGeo))
                {
                    // fractures don't exchange mass with barriers
                    if (scv.isOnFracture())
                        continue;

                    for (int pvIdx = 0; pvIdx < numEq; ++pvIdx)
                    {
                        if (isOnBoundary && bcTypes.isDirichlet(pvIdx))
                            continue;

                        auto evalFluxes = [&](Scalar priVar) {
                            elemSol[scv.localDofIndex()][pvIdx] = priVar;
                            deflectionHelper.deflect(elemSol, scv, this->problem());
                            return computeBarrierFlux(ip);
                        };

                        // derive the residuals numerically
                        static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
                        static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                        NumericDifferentiation::partialDerivative(
                            evalFluxes, elemSol[scv.localDofIndex()][pvIdx], partialDerivs, flux,
                            eps_(elemSol[scv.localDofIndex()][pvIdx], pvIdx), numDiffMethod
                        );

                        for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                            A[ip.insideScv.dofIndex()][scv.dofIndex()][eqIdx][pvIdx]
                                += partialDerivs[eqIdx];

                        deflectionHelper.restore(scv);
                        elemSol[scv.localDofIndex()][pvIdx] = curSol[scv.dofIndex()][pvIdx];
                    }
                }
            };

            auto insideElemSol = elementSolution(ip.insideFVGeometry.element(), ip.insideElemVolVars, ip.insideFVGeometry);
            auto outsideElemSol = elementSolution(ip.outsideFVGeometry.element(), ip.outsideElemVolVars, ip.outsideFVGeometry);
            addDerivatives(insideElemSol, ip.insideElemVolVars, ip.insideFVGeometry);
            addDerivatives(outsideElemSol, ip.outsideElemVolVars, ip.outsideFVGeometry);
        });
    }

    template<class SolutionVector, class Visitor>
    void visitBarriers_(const SolutionVector& curSol, const Visitor& visitor) const
    {
        auto insideFVGeometry = localView(this->gridGeometry());
        auto insideElemVolVars = localView(this->gridVariables().curGridVolVars());
        auto outsideFVGeometry = localView(this->gridGeometry());
        auto outsideElemVolVars = localView(this->gridVariables().curGridVolVars());
        for (const auto& e : elements(this->gridGeometry().gridView()))
        {
            insideFVGeometry.bindElement(e);
            insideElemVolVars.bindElement(e, insideFVGeometry, curSol);
            const auto& refElem = referenceElement(e);

            std::optional<typename GridGeometry::GridView::template Codim<0>::Entity::Geometry> elemGeo;
            const auto getElemGeo = [&] () {
                if (!elemGeo)
                    elemGeo = e.geometry();
                return elemGeo.value();
            };

            for (const auto& is : intersections(this->gridGeometry().gridView(), e))
                if (this->gridGeometry().isOnBarrier(e, is))
                {
                    const auto findNeighborScv = [&] (auto vIdx, const auto& fvGeo) -> const SubControlVolume& {
                        for (const auto& scv : scvs(fvGeo))
                            if (!scv.isOnFracture() &&
                                this->gridGeometry().vertexMapper().vertexIndex(fvGeo.element(), scv.localDofIndex(), dim) == vIdx)
                                return scv;
                        DUNE_THROW(Dune::InvalidStateException, "Could not find neighboring scv");
                    };

                    const auto& outside = is.outside();
                    const auto& normal = is.centerUnitOuterNormal();
                    outsideFVGeometry.bindElement(outside);
                    outsideElemVolVars.bindElement(outside, outsideFVGeometry, curSol);
                    visitIntegrationPoints_(getElemGeo(), is, [&] (int i, const auto& ip, const auto vol) {
                        const auto vIdx = this->gridGeometry().vertexMapper().vertexIndex(
                            e, refElem.subEntity(is.indexInInside(), 1, i, dim), dim
                        );
                        visitor(BarrierIntegrationPoint{
                            ip, normal, vol,
                            insideFVGeometry, outsideFVGeometry,
                            findNeighborScv(vIdx, insideFVGeometry), findNeighborScv(vIdx, outsideFVGeometry),
                            insideElemVolVars, outsideElemVolVars,
                        });
                    });
                }
        }
    }

    template<class ElemGeo, class Intersection, class Visitor>
    void visitIntegrationPoints_(const ElemGeo& elemGeo, const Intersection& is, const Visitor& visitor) const
    {
        for (int i = 0; i < is.geometry().corners(); ++i)
        {
            const auto [center, volume] = computeIntegrationPointData_(elemGeo, is.indexInInside(), i);
            visitor(i, center, volume);
        }
    }

    template<class Geometry>
    std::pair<GlobalPosition, typename GlobalPosition::value_type>
    computeIntegrationPointData_(const Geometry& geo, unsigned int facetIdx, int indexInFacet) const
    {
        using namespace Detail::Box;
        if constexpr (dim == 2)
        {
            const auto corners = subEntityKeyToCornerStorage<Dune::FieldVector<GlobalPosition, 2>>(
                geo, facetIdx, 1, ScvCorners<Dune::GeometryTypes::line>::keys[indexInFacet]
            );
            return std::make_pair(computeCenter_(corners), convexPolytopeVolume<dim-1>(Dune::GeometryTypes::line, [&] (auto i) {
                return corners[i];
            }));
        }
        else if constexpr (dim == 3)
        {
            if (geo.type() == Dune::GeometryTypes::tetrahedron)
            {
                const auto corners = subEntityKeyToCornerStorage<Dune::FieldVector<GlobalPosition, 4>>(
                    geo, facetIdx, 1, ScvCorners<Dune::GeometryTypes::triangle>::keys.at(indexInFacet)
                );
                return std::make_pair(computeCenter_(corners), convexPolytopeVolume<dim-1>(Dune::GeometryTypes::quadrilateral, [&] (auto i) {
                    return corners[i];
                }));
            }
            if (geo.type() == Dune::GeometryTypes::hexahedron)
            {
                const auto corners = subEntityKeyToCornerStorage<Dune::FieldVector<GlobalPosition, 4>>(
                    geo, facetIdx, 1, ScvCorners<Dune::GeometryTypes::quadrilateral>::keys[indexInFacet]
                );
                return std::make_pair(computeCenter_(corners), convexPolytopeVolume<dim-1>(Dune::GeometryTypes::quadrilateral, [&] (auto i) {
                    return corners[i];
                }));
            }
            else
                DUNE_THROW(Dune::NotImplemented, "Unsupported geometry type");
        }
    }

    template<class CornerStorage>
    auto computeCenter_(const CornerStorage& corners) const
    {
        std::decay_t<decltype(corners[0])> result(0.0);
        for (const auto& c : corners)
            result += c;
        result /= corners.size();
        return result;
    }
};

} // namespace Dumux

#endif
