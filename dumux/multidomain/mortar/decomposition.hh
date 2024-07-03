// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Class that represents a domain decomposition.
 *        Contains the subdomain & mortar grid geometries and connectivity information between them.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_DECOMPOSITION_HH
#define DUMUX_MULTIDOMAIN_MORTAR_DECOMPOSITION_HH

#include <vector>
#include <type_traits>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dumux/multidomain/glue.hh>

namespace Dumux::Mortar {

/*!
 * \file
 * \ingroup MultiDomain
 * \brief Class that represents a decomposition of a domain.
 * \tparam SubDomainGridGeometry The type used to represent subdomain grid geometries.
 ' \tparam MortarGridGeometries The type used to represent mortar grid geometries.
 */
template<typename SDGG, typename MGG>
class Decomposition
{
    using MDGlue = MultiDomainGlue<
        typename SDGG::GridView, typename MGG::GridView,
        typename SDGG::ElementMapper, typename MGG::ElementMapper
    >;

    struct Interface
    {
        std::size_t subDomainId;
        std::size_t mortarId;
        MDGlue glue;
    };

public:
    template<typename T>
    using Ptr = std::shared_ptr<T>;
    using SubDomainGridGeometry = SDGG;
    using MortarGridGeometry = MGG;
    using Glue = MDGlue;

    explicit Decomposition(std::vector<Ptr<SubDomainGridGeometry>>&& subDomainGGs,
                           std::vector<Ptr<MortarGridGeometry>>&& mortarGGs)
    : subDomainGridGeometries_{std::move(subDomainGGs)}
    , mortarGridGeometries_{std::move(mortarGGs)}
    {
        subDomainToInterfaceIds_.resize(subDomainGridGeometries_.size());
        mortarToInterfaceIds_.resize(mortarGridGeometries_.size());
        for (std::size_t sdId = 0; sdId < subDomainGridGeometries_.size(); ++sdId)
        {
            for (std::size_t mId = 0; mId < mortarGridGeometries_.size(); ++mId)
            {
                auto glue = makeGlue(subDomainGridGeometry(sdId), mortarGridGeometry(mId));

                if (glue.size() > 0)
                {
                    std::cout << "Found intersections between subdomain " << sdId << " and mortar " << mId << std::endl;
                    subDomainToInterfaceIds_[sdId].push_back(interfaces_.size());
                    mortarToInterfaceIds_[mId].push_back(interfaces_.size());
                    if (mortarToInterfaceIds_[mId].size() > 2)
                        DUNE_THROW(Dune::InvalidStateException, "Each mortar domain must only be between two subdomains");
                    interfaces_.emplace_back(Interface{sdId, mId, std::move(glue)});
                }
            }
        }
    }

    const SubDomainGridGeometry& subDomainGridGeometry(std::size_t id) const { return *subDomainGridGeometries_.at(id); }
    SubDomainGridGeometry& subDomainGridGeometry(std::size_t id) { return *subDomainGridGeometries_.at(id); }

    const MortarGridGeometry& mortarGridGeometry(std::size_t id) const { return *mortarGridGeometries_.at(id); }
    MortarGridGeometry& mortarGridGeometry(std::size_t id) { return *mortarGridGeometries_.at(id); }

private:
    std::vector<Ptr<SubDomainGridGeometry>> subDomainGridGeometries_;
    std::vector<Ptr<MortarGridGeometry>> mortarGridGeometries_;
    std::vector<Interface> interfaces_;

    std::vector<std::vector<std::size_t>> subDomainToInterfaceIds_;
    std::vector<std::vector<std::size_t>> mortarToInterfaceIds_;
};

} // end namespace Dumux::Mortar

#endif
