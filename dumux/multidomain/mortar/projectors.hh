// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief Projector implementations betwenn subdomains and mortars.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_PROJECTORS_HH
#define DUMUX_MULTIDOMAIN_MORTAR_PROJECTORS_HH

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dumux/discretization/projection/projector.hh>
#include "projectorinterface.hh"

namespace Dumux::Mortar {

/*!
 * \ingroup MultiDomain
 * \brief Projector based on an L2-projection.
 */
template<typename SolutionVector>
class L2Projector : public ProjectorInterface<SolutionVector>
{
    using Scalar = typename SolutionVector::field_type;
    using Projector = Dumux::Projector<Scalar>;

public:
    //! Export the type used for projection matrices
    using Matrix = typename Projector::Matrix;

    //! Constructor
    L2Projector(Matrix&& massMatrix, Matrix&& projMatrix)
    : projector_(std::move(massMatrix), std::move(projMatrix))
    {}

    //! Project the vector x into the target domain
    SolutionVector project(const SolutionVector& x) const override
    { return projector_.project(x); }

private:
    Projector projector_;
};

/*!
 * \ingroup MultiDomain
 * \brief Projector types supported by the factory.
 */
enum class ProjectorType { L2 };

/*!
 * \ingroup MultiDomain
 * \brief A factory for projector implementations.
 */
template<typename SolutionVector>
class ProjectorFactory
{
public:
    using ProjectorPtr = std::shared_ptr<ProjectorInterface<SolutionVector>>;

    struct ProjectorPair
    {
        ProjectorPtr toSubDomainTrace;
        ProjectorPtr toMortar;
    };

    ProjectorFactory(ProjectorType pt)
    : projectorType_{pt}
    {}

    /*!
     * \brief Create a forward/backwards projector pair.
     */
    template<typename SubDomainTraceBasis, typename MortarBasis, typename Glue>
    ProjectorPair make(const SubDomainTraceBasis& traceBasis,
                       const MortarBasis& mortarBasis,
                       const Glue& glue)
    {
        if (projectorType_ == ProjectorType::L2) {
            auto matrixPair = makeProjectionMatricesPair(traceBasis, mortarBasis, glue);
            return {
                .toSubDomainTrace = std::make_shared<L2Projector<SolutionVector>>(
                    std::move(matrixPair.second.first),
                    std::move(matrixPair.second.second)
                ),
                .toMortar = std::make_shared<L2Projector<SolutionVector>>(
                    std::move(matrixPair.first.first),
                    std::move(matrixPair.first.second)
                )
            };
        }
        DUNE_THROW(Dune::NotImplemented, "Unsupported projector type");
    }

private:
    ProjectorType projectorType_;
};

} // end namespace Dumux::Mortar

#endif
