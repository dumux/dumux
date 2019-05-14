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
#include <dumux/discretization/method.hh>
#include "common/projectorbase.hh"

namespace Dumux {

// namespace with some implementation details
namespace FluxProjectionDetail {

//! computes the interface pressure acting on a mortar element
//! This overload is the specialization for the box scheme used
template<class Scalar, DiscretizationMethod dm,
         class SubDomainGridGeometry, class SubDomainSolutionVector, class SubDomainIntegrationData,
         class MortarElement,
         std::enable_if_t<dm == DiscretizationMethod::box, int> = 0>
Scalar computePressureOnMortarElement(const SubDomainGridGeometry& gridGeometry,
                                      const SubDomainSolutionVector& x,
                                      const SubDomainIntegrationData& intData,
                                      const MortarElement& element)
{
    // integrate pressure acting from sub-domain on this element
    Scalar p = 0.0;
    for (const auto& dofToValuePair : intData.otherDofToWeightMap)
        p += x[dofToValuePair.first]*dofToValuePair.second;
    return p;
}

} // end namespace FluxProjectionDetail

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template< class MortarFEBasis, class MortarSolutionVector,
          class SubDomainGridGeometry, class SubDomainSolutionVector >
class MortarFluxProjector : public MortarProjectorBase<MortarFEBasis,
                                                       MortarSolutionVector,
                                                       SubDomainGridGeometry,
                                                       SubDomainSolutionVector>
{
    using ParentType = MortarProjectorBase< MortarFEBasis, MortarSolutionVector,
                                            SubDomainGridGeometry, SubDomainSolutionVector >;

    using Scalar = typename ParentType::Scalar;
    static constexpr DiscretizationMethod subDomainDM = SubDomainGridGeometry::discMethod;
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

            p[eIdx] = FluxProjectionDetail:: computePressureOnMortarElement<Scalar, subDomainDM>(*this->subDomainGridGeometry_,
                                                                                                 *this->subDomainSolution_,
                                                                                                 it->second,
                                                                                                 element);
        }

        return p;
    }
};

} // end namespace Dumux

#endif
