// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief The specialized Dumux tag for the ISTL registry including only ISTL solvers compatible
 *        with MultiTypeBlockMatrix
 */
#ifndef DUMUX_LINEAR_ISTL_SOLVERS_MULTITYPE_HH
#define DUMUX_LINEAR_ISTL_SOLVERS_MULTITYPE_HH

#include <dune/common/version.hh>

#include <dune/istl/solvers.hh>
#include <dune/istl/cholmod.hh>
#include <dune/istl/umfpack.hh>

#include <dumux/linear/istlsolverregistry.hh>

namespace Dumux {

DUMUX_REGISTER_SOLVER("loopsolver", Dumux::MultiTypeBlockMatrixSolverTag, defaultIterativeSolverCreator<Dune::LoopSolver>());
DUMUX_REGISTER_SOLVER("gradientsolver", Dumux::MultiTypeBlockMatrixSolverTag, defaultIterativeSolverCreator<Dune::GradientSolver>());
DUMUX_REGISTER_SOLVER("cgsolver", Dumux::MultiTypeBlockMatrixSolverTag, defaultIterativeSolverCreator<Dune::CGSolver>());
DUMUX_REGISTER_SOLVER("bicgstabsolver", Dumux::MultiTypeBlockMatrixSolverTag, defaultIterativeSolverCreator<Dune::BiCGSTABSolver>());
DUMUX_REGISTER_SOLVER("minressolver", Dumux::MultiTypeBlockMatrixSolverTag, defaultIterativeSolverCreator<Dune::MINRESSolver>());
DUMUX_REGISTER_SOLVER("restartedgmressolver", Dumux::MultiTypeBlockMatrixSolverTag, defaultIterativeSolverCreator<Dune::RestartedGMResSolver>());
DUMUX_REGISTER_SOLVER("restartedflexiblegmressolver", Dumux::MultiTypeBlockMatrixSolverTag, defaultIterativeSolverCreator<Dune::RestartedFlexibleGMResSolver>());
DUMUX_REGISTER_SOLVER("generalizedpcgsolver", Dumux::MultiTypeBlockMatrixSolverTag, defaultIterativeSolverCreator<Dune::GeneralizedPCGSolver>());
DUMUX_REGISTER_SOLVER("restartedfcgsolver", Dumux::MultiTypeBlockMatrixSolverTag, defaultIterativeSolverCreator<Dune::RestartedFCGSolver>());
DUMUX_REGISTER_SOLVER("completefcgsolver", Dumux::MultiTypeBlockMatrixSolverTag, defaultIterativeSolverCreator<Dune::CompleteFCGSolver>());
#if HAVE_SUITESPARSE_CHOLMOD
DUMUX_REGISTER_SOLVER("cholmod", Dumux::MultiTypeBlockMatrixSolverTag,
                      [](auto opTraits, const auto& op, const Dune::ParameterTree& config)
                      -> std::shared_ptr<typename decltype(opTraits)::solver_type>
                      {
                        using OpTraits = decltype(opTraits);
                        using M = typename OpTraits::matrix_type;
                        using D = typename OpTraits::domain_type;
                        // works only for sequential operators
                        if constexpr (OpTraits::isParallel){
                          if(opTraits.getCommOrThrow(op).communicator().size() > 1)
                            DUNE_THROW(Dune::InvalidStateException, "CholMod works only for sequential operators.");
                        }
                        if constexpr (OpTraits::isAssembled &&
                                      // check whether the Matrix field_type is double or float
                                      (std::is_same_v<typename FieldTraits<D>::field_type, double> ||
                                      std::is_same_v<typename FieldTraits<D>::field_type, float>)){
                          const auto& A = opTraits.getAssembledOpOrThrow(op);
                          const M& mat = A->getmat();
                          auto solver = std::make_shared<Dune::Cholmod<D>>();
                          solver->setMatrix(mat);
                          return solver;
                        }
                        DUNE_THROW(UnsupportedType,
                                   "Unsupported Type in Cholmod (only double and float supported)");
                      });
#endif // HAVE_SUITESPARSE_CHOLMOD
#if HAVE_SUITESPARSE_UMFPACK
DUMUX_REGISTER_SOLVER("umfpack", Dumux::MultiTypeBlockMatrixSolverTag,
                      [](auto opTraits, const auto& op, const Dune::ParameterTree& config)
                      -> std::shared_ptr<typename decltype(opTraits)::solver_type>
                      {
                        using OpTraits = decltype(opTraits);
                        // works only for sequential operators
                        if constexpr (OpTraits::isParallel){
                          if(opTraits.getCommOrThrow(op).communicator().size() > 1)
                            DUNE_THROW(Dune::InvalidStateException, "UMFPack works only for sequential operators.");
                        }
                        if constexpr (OpTraits::isAssembled){
                          using M = typename OpTraits::matrix_type;
                          // check if UMFPack<M>* is convertible to
                          // InverseOperator*. This checks compatibility of the
                          // domain and range types
                          if constexpr (UMFPackImpl::isValidBlock<OpTraits>::value) {
                            const auto& A = opTraits.getAssembledOpOrThrow(op);
                            const M& mat = A->getmat();
                            int verbose = config.get("verbose", 0);
                            return std::make_shared<Dune::UMFPack<M>>(mat,verbose);
                          }
                        }
                        DUNE_THROW(UnsupportedType,
                                   "Unsupported Type in UMFPack (only double and std::complex<double> supported)");
                        return nullptr;
                      });
#endif // HAVE_SUITESPARSE_UMFPACK
} // end namespace Dumux

#endif // DUMUX_LINEAR_ISTL_SOLVERS_MULTITYPE_HH
