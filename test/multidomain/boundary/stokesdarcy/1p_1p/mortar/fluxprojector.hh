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
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
#ifndef DUMUX_MORTAR_FLUX_PROJECTOR_HH
#define DUMUX_MORTAR_FLUX_PROJECTOR_HH

#include <string>
#include <unordered_map>

#include <dune/common/exceptions.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include "common/projectorbase.hh"

namespace Dumux {

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template< class MortarFEBasis, class MortarSolutionVector,
          class SubDomainGridGeometry, class SubDomainSolutionVector,
          DiscretizationMethod subDomainDM = SubDomainGridGeometry::discMethod>
class MortarFluxProjector;

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template< class MortarFEBasis, class MortarSolutionVector,
       class SubDomainGridGeometry, class SubDomainSolutionVector >
class MortarFluxProjector< MortarFEBasis, MortarSolutionVector,
                           SubDomainGridGeometry, SubDomainSolutionVector,
                           DiscretizationMethod::box>
: public MortarProjectorBase< MortarFEBasis, MortarSolutionVector,
                              SubDomainGridGeometry, SubDomainSolutionVector >
{
    using ParentType = MortarProjectorBase< MortarFEBasis, MortarSolutionVector,
                                            SubDomainGridGeometry, SubDomainSolutionVector >;

    using Scalar = typename ParentType::Scalar;
    static constexpr DiscretizationMethod subDomainDM = DiscretizationMethod::box;
public:

    //! The constructor
    MortarFluxProjector(std::shared_ptr<const MortarFEBasis> mortarFEBasis,
                        std::shared_ptr<const SubDomainGridGeometry> subDomainGridGeometry,
                        const std::string& paramGroup = "")
    : ParentType(mortarFEBasis, subDomainGridGeometry, paramGroup)
    {}

    //! projects the sub-domain interface pressures to mortar space
    MortarSolutionVector projectInterfacePressures() const
    {
        MortarSolutionVector p;
        p.resize(this->mortarFEBasis_->gridView().size(0));

        for (const auto& element : elements(this->mortarFEBasis_->gridView()))
        {
            const auto eIdx = this->mortarElementMapper_.index(element);
            auto it = this->subDomainToMortarProjection_.find(eIdx);

            // only proceed if there is an entry
            if (it == this->subDomainToMortarProjection_.end())
                continue;

            Scalar elemPressure = 0.0;
            for (const auto& dofToValuePair : it->second.otherDofToWeightMap)
                elemPressure += (*this->subDomainSolution_)[dofToValuePair.first]*dofToValuePair.second;
            p[eIdx] = elemPressure;
        }

        return p;
    }
};

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template< class MortarFEBasis, class MortarSolutionVector,
          class SubDomainGridGeometry, class SubDomainSolutionVector >
class MortarFluxProjector< MortarFEBasis, MortarSolutionVector,
                           SubDomainGridGeometry, SubDomainSolutionVector,
                           DiscretizationMethod::cctpfa>
: public MortarProjectorBase< MortarFEBasis, MortarSolutionVector,
                              SubDomainGridGeometry, SubDomainSolutionVector >
{
    using ParentType = MortarProjectorBase< MortarFEBasis, MortarSolutionVector,
                                            SubDomainGridGeometry, SubDomainSolutionVector >;

    using SubDomainGridView = typename SubDomainGridGeometry::GridView;
    using SubDomainGridIndex = typename IndexTraits<SubDomainGridView>::GridIndex;
    static constexpr DiscretizationMethod subDomainDM = DiscretizationMethod::box;
public:

    //! export type used for scalar values
    using typename ParentType::Scalar;

    //! The constructor
    MortarFluxProjector(std::shared_ptr<const MortarFEBasis> mortarFEBasis,
                        std::shared_ptr<const SubDomainGridGeometry> subDomainGridGeometry,
                        const std::string& paramGroup = "")
    : ParentType(mortarFEBasis, subDomainGridGeometry, paramGroup)
    {
        // save to each sub-domain element the scvf index that couples to mortar
        for (const auto& mortarElement : elements(mortarFEBasis->gridView()))
        {
            const auto mortarElementIdx = this->mortarElementMapper_.index(mortarElement);
            auto it = this->subDomainToMortarProjection_.find(mortarElementIdx);
            if (it == this->subDomainToMortarProjection_.end())
                continue;

            const auto mortarEG = mortarElement.geometry();
            const auto mortarElementCenter = mortarEG.center();

            auto mortarElementEdge = mortarEG.corner(1) - mortarEG.corner(0);
            mortarElementEdge /= mortarElementEdge.two_norm();

            const auto& coupledSDElementIndices = it->second.coupledElementIndices;
            for (const auto& sdElementIdx : coupledSDElementIndices)
            {
                auto& mapEntry = subDomainElementToCoupledScvfIndex_[sdElementIdx];
                const auto subDomainElement = subDomainGridGeometry.element(sdElementIdx);

                auto fvGeometry = localView(subDomainGridGeometry);
                fvGeometry.bind(subDomainElement);

                using std::abs;
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    // mortar edge is normal to scvf it is a candidate
                    if ( abs(scvf.unitOuterNormal()*mortarElementEdge) > 1e-8 )
                        continue;

                    auto d = scvf.center() - mortarElementCenter;
                    d /= d.two_norm();

                    // it is a coupled face
                    if ( abs(scvf.unitOuterNormal()*d) < 1e-8 )
                    {
                        if (!scvf.boundary())
                            DUNE_THROW(Dune::InvalidStateException, "Found non-boundary scvf to be coupling scvf!");
                        else if (mapEntry.size() == 0)
                            mapEntry.push_back(scvf.index());
                        else if (mapEntry[0] != scvf.index())
                            DUNE_THROW(Dune::InvalidStateException, "Found two coupling scvfs in sub-domain!");
                        else
                            DUNE_THROW(Dune::InvalidStateException, "Unexpected code behaviour");
                    }
                }
            }
        }
    }

    //! projects the sub-domain interface pressures to mortar space
    MortarSolutionVector projectInterfacePressures() const
    {
        DUNE_THROW(Dune::NotImplemented, "Tpfa reconstruction not yet implemented");

        // first compute all interface pressures in sub-domain
        std::unordered_map<SubDomainGridIndex, Scalar> sdPressures;

        for (const auto& entry : subDomainElementToCoupledScvfIndex_)
        {
            const auto sdElemIdx = entry.first;
            const auto sdScvfIdx = entry.second;
            const auto sdElement = this->subDomainGridGeometry_->element(sdElemIdx);

            // compute mortar flux for this element
            const auto flux = this->integrateMortarVariable(sdElement);

            auto fvGeometry = localView(*this->subDomainGridGeometry_);
            fvGeometry.bindElement(sdElement);

            const auto& scvf = fvGeometry.scvf(sdScvfIdx);
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto elemSol = elementSolution(sdElement,
                                                 *this->subDomainSolution_,
                                                 *this->subDomainGridGeometry_);

            // TODO: how to get K/rho/mu???
            // const auto exFactor = problem.extrusionFactor(element, insideScv, elemSol);
            // const auto k = problem.spatialParams().permeability(element, insideScv, elemSol);
            // const auto ti = computeTpfaTransmissibility(scvf, insideScv, k, exFactor);
        }
    }

private:
    std::unordered_map< SubDomainGridIndex, std::vector<SubDomainGridIndex> > subDomainElementToCoupledScvfIndex_;
};

} // end namespace Dumux

#endif
