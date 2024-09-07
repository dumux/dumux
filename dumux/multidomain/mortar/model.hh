// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TODO
 * \brief TODO
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_MODEL_HH
#define DUMUX_MULTIDOMAIN_MORTAR_MODEL_HH

#include <vector>
#include <memory>
#include <limits>
#include <concepts>

#include <dumux/multidomain/glue.hh>

#include "projectorinterface.hh"
#include "trace.hh"

namespace Dumux::Mortar {

template<typename SolutionVector>
struct SubDomain
{
    virtual void solve() = 0;
    virtual void setMortarTrace(std::size_t, SolutionVector&&) = 0;
    virtual SolutionVector assembleTraceWith(std::size_t) const = 0;
};

template<typename SolutionVector>
class Projector
{
    using Interface = ProjectorInterface<SolutionVector>;
public:

    template<std::derived_from<Interface> ToTrace,
             std::derived_from<Interface> FromTrace>
    Projector(ToTrace&& toTrace, FromTrace&& fromTrace)
    : toSubDomainTrace_{std::move(toTrace)}
    , fromSubDomainTrace_{std::move(fromTrace)}
    {}

    SolutionVector toSubDomainTrace(const SolutionVector& x) const
    { return toSubDomainTrace_->project(x); }

    SolutionVector fromSubDomainTrace(const SolutionVector& x) const
    { return fromSubDomainTrace_->project(x); }

private:
    std::unique_ptr<Interface> toSubDomainTrace_;
    std::unique_ptr<Interface> fromSubDomainTrace_;
};

template<typename SolutionVector>
struct Mortar
{
    virtual const SolutionVector& solution() const = 0;
};

struct Interface
{
    struct Side
    {
        std::size_t subDomainIndex;
        std::size_t projectorIndex;
    };

    std::size_t mortarIndex;
    Side positive;
    Side negative;
};

template<typename SolutionVector>
class Model
{
 public:
    void splitMortarSolution(const SolutionVector& x, std::vector<SolutionVector>& result) const
    {
        for (std::size_t id = 0; id < mortars_.size(); ++id)
        {
            const auto& x = mortars_[id]->solution();
            result[id].resize(x.size());
            std::ranges::copy(x, result[id].begin());
        }
    }

    void projectMortarToSubDomains(const std::vector<SolutionVector>& mortarSolutions) const
    {
        for (const auto& interface : interfaces_)
        {
            subDomains_[interface.positive.subDomainIndex].setMortarTrace(
                interface.mortarIndex,
                projectors_[interface.positive.projectorIndex].toSubDomainTrace(mortarSolutions[interface.mortarIndex])
            );
            subDomains_[interface.negative.subDomainIndex].setMortarTrace(
                interface.mortarIndex,
                projectors_[interface.negative.projectorIndex].toSubDomainTrace(mortarSolutions[interface.mortarIndex])
            );
        }
    }

    void solveSubDomains() const
    {
        for (const auto& sdPtr : subDomains_)
            sdPtr->solve();
    }

    void computeMortarJump(std::vector<SolutionVector>& result) const
    {
        for (const auto& interface : interfaces_)
        {
            const auto& positive = subDomains_[interface.positive.subDomainIndex];
            const auto& negative = subDomains_[interface.negative.subDomainIndex];
            const auto posTrace = positive.assembleTraceWith(interface.mortarIndex);
            const auto negTrace = negative.assembleTraceWith(interface.mortarIndex);
            result[interface.mortarIndex] = projectors_[interface.positive.projectorIndex].fromSubDomainTrace(posTrace);
            result[interface.mortarIndex] -= projectors_[interface.negative.projectorIndex].fromSubDomainTrace(negTrace);
        }
    }

    void mergeMortarSolutions(const std::vector<SolutionVector>& in, SolutionVector& result) const
    {
        result.resize(
            std::accumulate(
                mortars_.begin(), mortars_.end(), std::size_t{1},
                [] (std::size_t current, const auto& mortar) {
                    return current + mortar->solution().size();
                }
            )
        );

        std::size_t offset = 0;
        for (std::size_t id = 0; id < mortars_.size(); ++id)
        {
            const auto& x = mortars_[id]->solution();
            std::ranges::copy(x, result.begin() + offset);
            offset += x.size();
        }
    }

 private:
    std::vector<Interface> interfaces_;
    std::vector<std::unique_ptr<Projector<SolutionVector>>> projectors_;
    std::vector<std::unique_ptr<SubDomain<SolutionVector>>> subDomains_;
    std::vector<std::unique_ptr<Mortar<SolutionVector>>> mortars_;
};

// TODO: Distinguish mortar and SD solution vector
// TODO: Template on MortarGridView and hardcode FEGridGeometry?
template<typename SolutionVector, typename MortarGridGeometry>
class ModelFactory
{
    static constexpr std::size_t invalidId = std::numeric_limits<std::size_t>::max();

    class MortarImpl : public Mortar<SolutionVector> {
     public:
        MortarImpl(std::shared_ptr<MortarGridGeometry> gg)
        : gridGeometry_{std::move(gg)}
        , solution_{}
        {
            solution_.resize(gridGeometry_->numDofs());
            solution_ = 0;
        }

        const MortarGridGeometry& gridGeometry() const
        { return *gridGeometry_; }

        const SolutionVector& solution() const override
        { return solution_; }

     private:
        std::shared_ptr<MortarGridGeometry> gridGeometry_;
        SolutionVector solution_;
    };

    // type-erased subdomain implementation
    template<typename SD>
    class SubDomainImpl : public SubDomain<SolutionVector> {
     public:
        SubDomainImpl(std::shared_ptr<SD> sd)
        : subDomain_{std::move(sd)}
        {}

        void solve() override
        { subDomain_->solve(); }

        void setMortarTrace(std::size_t i, SolutionVector&& x) override
        { subDomain_->setMortarTrace(i, std::move(x)); }

        SolutionVector assembleTraceWith(std::size_t i) const override
        { return subDomain_->assembleTraceWith(i); }

     private:
        std::shared_ptr<SD> subDomain_;
    };

 public:
    void insertMortar(std::shared_ptr<MortarGridGeometry> gridGeometry) {
        if (!subDomains_.empty())
            DUNE_THROW(Dune::InvalidStateException, "Mortars have to be inserted prior to the subdomains");
        mortars_.emplace_back(MortarImpl{std::move(gridGeometry)});
        mortarToInterface_.push_back(invalidId);
    }

    template<typename SD, typename InterfaceCallBack>  // TODO: SD concept?; TODO: Callback concept?
    void insertSubDomain(std::shared_ptr<SD> subDomain, InterfaceCallBack&& callback) {
        FVTrace trace{subDomain->gridVariables(), [] (auto&&...) { return true; }};

        using TraceGridView = std::remove_cvref_t<decltype(trace.gridView())>;
        using TraceElementMapper = std::remove_cvref_t<decltype(trace.elementMapper())>;
        using TraceEntitySet = GridViewGeometricEntitySet<TraceGridView, 0, TraceElementMapper>;
        BoundingBoxTree<TraceEntitySet> traceBBoxTree{std::make_shared<TraceEntitySet>(trace.gridView())};

        subDomainToInterfaces_.push_back({});
        for (std::size_t mortarId = 0; mortarId < mortars_.size(); ++mortarId)
        {
            MultiDomainGlue<
                TraceGridView, typename MortarGridGeometry::GridView,
                TraceElementMapper, typename MortarGridGeometry::ElementMapper
            > glue;

            glue.build(traceBBoxTree, mortars_[mortarId].gridGeometry().boundingBoxTree());
            if (glue.size() > 0)
            {
                std::cout << "Subdomain " << subDomains_.size() << "  -> mortar " << mortarId << std::endl;

                FVTrace subTrace{subDomain->gridVariables(), [&] (const auto& is) {
                    return !intersectingEntities(is.geometry(), mortars_[mortarId].gridGeometry().boundingBoxTree()).empty();
                }};

                callback(subDomain, std::move(subTrace), mortarId);
                registerMapping_(subDomains_.size(), mortarId);
            }
        }
        subDomains_.emplace_back(std::make_unique<SubDomainImpl<SD>>(std::move(subDomain)));
    }

 private:
    void registerMapping_(std::size_t subDomainIndex, std::size_t mortarIndex)
    {
        const auto interfaceId = mortarToInterface_.at(mortarIndex);
        if (interfaceId != invalidId)
        {
            auto& interface = interfaces_.at(interfaceId);
            if (interface.negative.subDomainIndex != invalidId)
                DUNE_THROW(Dune::InvalidStateException, "Negative side already registered for this mortar");
            interface.negative.subDomainIndex = subDomainIndex;
        }
        else
        {
            mortarToInterface_.at(mortarIndex) = interfaces_.size();
            interfaces_.emplace_back(Interface{
                .mortarIndex = mortarIndex,
                .positive = {.subDomainIndex = subDomainIndex, .projectorIndex = invalidId},
                .negative = {invalidId, invalidId}
            });
        }
    }

    std::vector<Interface> interfaces_;
    std::vector<std::unique_ptr<Projector<SolutionVector>>> projectors_;
    std::vector<std::unique_ptr<SubDomain<SolutionVector>>> subDomains_;
    std::vector<MortarImpl> mortars_;
    std::vector<std::vector<std::size_t>> subDomainToInterfaces_;
    std::vector<std::size_t> mortarToInterface_;
};

} // end namespace Dumux::Mortar

#endif
