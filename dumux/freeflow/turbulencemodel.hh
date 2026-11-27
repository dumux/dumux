// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \brief The available free flow turbulence models in Dumux
 */
#ifndef DUMUX_FREEFLOW_TURBULENCE_MODEL_HH
#define DUMUX_FREEFLOW_TURBULENCE_MODEL_HH

#include <string>
#include <dumux/common/exceptions.hh>

namespace Dumux {

/*!
 * \brief The available free flow turbulence models in Dumux
 * \ingroup FreeflowModels
 * \note Use none for plain (Navier-) Stokes models (DNS)
 */
enum class TurbulenceModel
{
    none, zeroeq, oneeq, kepsilon, lowrekepsilon, komega, sst
};

/*!
 * \brief The number of extra transport equations a given turbulence model adds
 *        on top of the plain (Navier-)Stokes mass and momentum balance equations.
 */
constexpr unsigned int numTurbulenceEq(TurbulenceModel model)
{
    if (model == TurbulenceModel::none || model == TurbulenceModel::zeroeq)
        return 0;
    else if (model == TurbulenceModel::oneeq)
        return 1;
    else
        return 2;
}

/*!
 * \brief return the name of the turbulence model
 */
inline std::string turbulenceModelToString(TurbulenceModel turbulenceModel)
{
    switch (turbulenceModel)
    {
        case TurbulenceModel::none: return "No_TurbModel";
        case TurbulenceModel::zeroeq: return "ZeroEq";
        case TurbulenceModel::oneeq: return "OneEq";
        case TurbulenceModel::kepsilon: return "KEpsilon";
        case TurbulenceModel::lowrekepsilon: return "LowReKEpsilon";
        case TurbulenceModel::komega: return "KOmega";
        case TurbulenceModel::sst: return "KOmegaSST";
        default: return "Invalid"; // should never be reached
    }
}

/*!
 * \brief The available variations of the SST turbulence model
 * \ingroup FreeflowModels
 */
enum class SSTModel
{ BSL, SST };

/*!
 * \brief return the name of the SST model variant as a string
 */
inline std::string sstModelToString(SSTModel sstModel)
{
    switch (sstModel)
    {
        case SSTModel::BSL: return "BSL";
        case SSTModel::SST: return "SST";
        default: return "Invalid";
    }
}

/*!
 * \brief Convenience function to convert user input given as std::string
 *        to the corresponding enum class used for choosing the SST model variant
 */
inline SSTModel sstModelFromString(const std::string& sstModel)
{
    if (sstModel == "BSL") return SSTModel::BSL;
    if (sstModel == "SST") return SSTModel::SST;
    DUNE_THROW(ParameterException, "\nThis SST model approach: \"" << sstModel << "\" is not implemented.\n"
                                   << "The available SST models are as follows: \n"
                                   << sstModelToString(SSTModel::BSL) << ": The baseline SST model\n"
                                   << sstModelToString(SSTModel::SST) << ": The full standard SST model");
}

} // end namespace Dumux

#endif
