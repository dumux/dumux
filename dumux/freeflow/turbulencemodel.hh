// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup FreeflowModels
 * \brief The available free flow turbulence models in Dumux
 */
#ifndef DUMUX_FREEFLOW_TURBLENCEMODEL_HH
#define DUMUX_FREEFLOW_TURBLENCEMODEL_HH

#include <string>

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

    constexpr unsigned int numTurbulenceEq(TurbulenceModel model)
    {
        if (model == TurbulenceModel::none || model == TurbulenceModel::zeroeq)
            return 0;
        else if (model == TurbulenceModel::oneeq)
            return 1;
        else
            return 2;
    }

    /**
     * \brief return the name of the Turbulence Model
     */
    std::string turbulenceModelToString(TurbulenceModel turbulenceModel)
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
     * \brief The available variations of the SST Turbulence Model
     * \ingroup SSTModel
     */
    enum class SSTModel
    { BSL, SST };

    /**
     * \brief return the name of the sst Model as a string
     */
    std::string sstModelToString(SSTModel sstModel)
    {
        switch (sstModel)
        {
            case SSTModel::BSL: return "BSL";
            case SSTModel::SST: return "SST";
            default: return "Invalid";
        }
    }

    /**
     * \brief Convenience function to convert user input given as std::string
     *        to the corresponding enum class used for choosing the SST Model
     */
    SSTModel sstModelFromString(const std::string& sstModel)
    {
        if (sstModel == "BSL") return SSTModel::BSL;
        if (sstModel == "SST") return SSTModel::SST;
        DUNE_THROW(ParameterException, "\nThis SST Model approach : \"" << sstModel << "\" is not implemented.\n"
                                       << "The available SST models are as follows: \n"
                                       << sstModelToString(SSTModel::BSL) << ": The Baseline SST Model n\n"
                                       << sstModelToString(SSTModel::SST) << ": The full standard SST Model");
    }

} // end namespace Dumux

#endif
