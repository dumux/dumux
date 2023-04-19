// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Index helpers for the free-flow/porous-medium-flow coupling.
 */

#ifndef DUMUX_MD_FREEFLOW_POROUSMEDIUM_INDEX_HELPER_HH
#define DUMUX_MD_FREEFLOW_POROUSMEDIUM_INDEX_HELPER_HH

#include <dune/common/indices.hh>

namespace Dumux::FreeFlowPorousMediumCoupling {

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model.
 * \tparam freeFlowIdx The domain index of the free-flow model.
 * \tparam porousMediumIndex The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 * \tparam hasAdapter Specifies whether an adapter class for the fluidsystem is used.
 */
template<std::size_t freeFlowIdx, std::size_t porousMediumIndex, class FFFS, bool hasAdapter>
struct IndexHelper;

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model. Specialization for the case that no adapter is used.
 * \tparam freeFlowIdx The domain index of the free-flow model.
 * \tparam porousMediumIndex The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 */
template<std::size_t freeFlowIdx, std::size_t porousMediumIndex, class FFFS>
struct IndexHelper<freeFlowIdx, porousMediumIndex, FFFS, false>
{
    /*!
     * \brief No adapter is used, just return the input index.
     */
    template<std::size_t i>
    static constexpr auto couplingPhaseIdx(Dune::index_constant<i>, int coupledPhaseIdx = 0)
    { return coupledPhaseIdx; }

    /*!
     * \brief No adapter is used, just return the input index.
     */
    template<std::size_t i>
    static constexpr auto couplingCompIdx(Dune::index_constant<i>, int coupledCompdIdx)
    { return coupledCompdIdx; }
};

/*!
 * \ingroup FreeFlowPorousMediumCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model. Specialization for the case that a adapter is used.
 * \tparam freeFlowIdx The domain index of the free-flow model.
 * \tparam porousMediumIndex The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 */
template<std::size_t freeFlowIdx, std::size_t porousMediumIndex, class FFFS>
struct IndexHelper<freeFlowIdx, porousMediumIndex, FFFS, true>
{
    /*!
     * \brief The free-flow model always uses phase index 0.
     */
    static constexpr int couplingPhaseIdx(Dune::index_constant<freeFlowIdx>, int coupledPhaseIdx = 0)
    { return 0; }

    /*!
     * \brief The phase index of the porous-medium-flow model is given by the adapter fluidsystem (i.e., user input).
     */
    static constexpr auto couplingPhaseIdx(Dune::index_constant<porousMediumIndex>, int coupledPhaseIdx = 0)
    { return FFFS::multiphaseFluidsystemPhaseIdx; }

    /*!
     * \brief The free-flow model does not need any change of the component index.
     */
    static constexpr auto couplingCompIdx(Dune::index_constant<freeFlowIdx>, int coupledCompIdx)
    { return coupledCompIdx; }

    /*!
     * \brief The component index of the porous-medium-flow model is mapped by the adapter fluidsystem.
     */
    static constexpr auto couplingCompIdx(Dune::index_constant<porousMediumIndex>, int coupledCompIdx)
    { return FFFS::compIdx(coupledCompIdx); }
};

} // end namespace Dumux::FreeFlowPorousMediumCoupling

#endif
