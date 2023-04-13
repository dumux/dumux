// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup Nonlinear
 * \brief This class provides the infrastructure to write the
 *        convergence behaviour of the Newton method for
 *        multidomain simulations into a VTK file.
 */
#ifndef DUMUX_MULTIDOMAIN_NEWTON_CONVERGENCE_WRITER_HH
#define DUMUX_MULTIDOMAIN_NEWTON_CONVERGENCE_WRITER_HH

#include <string>
#include <dune/common/hybridutilities.hh>
#include <dumux/nonlinear/newtonconvergencewriter.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup Nonlinear
 * \brief Writes the intermediate solutions for every Newton iteration
 * \note To use this create a shared_ptr to an instance of this class in the main file
 *       and pass it to newton.solve(x, convergencewriter). You can use the reset method
 *       to write out multiple Newton solves with a unique id, if you don't call use all
 *       Newton iterations just come after each other in the pvd file.
 */
template <class MDTraits>
class MultiDomainNewtonConvergenceWriter : public ConvergenceWriterInterface<typename MDTraits::SolutionVector,typename MDTraits::ResidualVector>
{
    template<std::size_t id>
    using SubDomainGridGeometry = typename MDTraits::template SubDomain<id>::GridGeometry;

    using GridGeometryTuple = typename MDTraits::template TupleOfSharedPtrConst<SubDomainGridGeometry>;

    using SolutionVector = typename MDTraits::SolutionVector;
    using ResidualVector = typename MDTraits::ResidualVector;

    template<std::size_t id>
    using SubDomainSolutionVector = typename MDTraits::template SubDomain<id>::SolutionVector;
    template<std::size_t id>
    using SubDomainResidualVector = typename MDTraits::template SubDomain<id>::ResidualVector;

    template<std::size_t id>
    using SubDomainNewtonConvergenceWriter = NewtonConvergenceWriter<SubDomainGridGeometry<id>, SubDomainSolutionVector<id>, SubDomainResidualVector<id>>;

    using ConvergenceWriterPtrTuple = typename MDTraits::template TupleOfSharedPtr<SubDomainNewtonConvergenceWriter>;

public:
    /*!
     * \brief Constructor
     * \param gridGeometryTuple A tuple of grid geometries
     * \param name Base name of the vtk output
     */
    MultiDomainNewtonConvergenceWriter(GridGeometryTuple gridGeometryTuple,
                                       const std::string& name = "newton_convergence")
    : gridGeometryTuple_(std::move(gridGeometryTuple))
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<MDTraits::numSubDomains>{}, [&](auto&& id)
        {
            using ConvWriter = SubDomainNewtonConvergenceWriter<std::decay_t<decltype(id)>::value>;
            elementAt(convergenceWriterPtrTuple_, id) = std::make_shared<ConvWriter>(*elementAt(gridGeometryTuple, id), name + "_domain_" + std::to_string(id));
        });
    }

    //! Resizes the output fields. This has to be called whenever the grid changes
    void resize()
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<MDTraits::numSubDomains>{}, [&](auto&& id)
        {
            elementAt(convergenceWriterPtrTuple_, id)->resize();
        });
    }

    //! Reset the convergence writer for a possible next Newton step
    //! You may set a different id in case you don't want the output to be overwritten by the next step
    void reset(std::size_t newId = 0UL)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<MDTraits::numSubDomains>{}, [&](auto&& id)
        {
            elementAt(convergenceWriterPtrTuple_, id)->reset(newId);
        });
    }

    void write(const SolutionVector& uLastIter,
               const ResidualVector& deltaU,
               const ResidualVector& residual) override
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<MDTraits::numSubDomains>{}, [&](auto&& id)
        {
            constexpr auto i = std::decay_t<decltype(id)>{};
            elementAt(convergenceWriterPtrTuple_, id)->write(uLastIter[i], deltaU[i], residual[i]);
        });
    }

private:
    GridGeometryTuple gridGeometryTuple_;
    ConvergenceWriterPtrTuple convergenceWriterPtrTuple_;
};

} // end namespace Dumux

#endif
