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
 * \ingroup Discretization
 * \brief Class that contains the element-local (or element stencil-local) data
 *        required to evaluate the terms of discrete equations.
 */
#ifndef DUMUX_LOCAL_CONTEXT_HH
#define DUMUX_LOCAL_CONTEXT_HH

#include <optional>
#include <dune/common/exceptions.hh>
#include <dumux/timestepping/timelevel.hh>

namespace Dumux {

class EmptyCouplingContext {};

/*!
 * \ingroup Discretization
 * \brief TODO: Doc me
 */
template<class ElementSolution,
         class ElementVariables,
         class CouplingContext = EmptyCouplingContext>
class LocalContext
{
    using Scalar = typename ElementSolution::PrimaryVariables::value_type;
    using TimeLevel = Dumux::TimeLevel<Scalar>;

public:
    //! TODO: Doc me
    void setElementSolution(const ElementSolution& elemSol)
    { elemSol_ = &elemSol; }

    //! TODO: Doc me
    void setElementVariables(const ElementVariables& elemVars)
    { elemVars_ = &elemVars; }

    //! TODO: Doc me
    void setTimeLevel(const TimeLevel& timeLevel)
    { timeLevel_.emplace(timeLevel); }

    //! TODO: Doc me
    void setCouplingContext(const CouplingContext& couplingContext)
    { couplingContext_ = &couplingContext; }

    //! TODO: Doc me
    const ElementSolution& elementSolution() const
    {
        if (!elemSol_)
            DUNE_THROW(Dune::InvalidStateException, "Element solution not available!");
        return *elemSol_;
    }

    //! TODO: Doc me
    const ElementVariables& elementVariables() const
    {
        if (!elemVars_)
            DUNE_THROW(Dune::InvalidStateException, "Element variables not available!");
        return *elemVars_;
    }

    //! TODO: Doc me
    const TimeLevel& timeLevel() const
    {
        if (!timeLevel_)
            DUNE_THROW(Dune::InvalidStateException, "Time level not available!");
        return *timeLevel_;
    }

    //! TODO: Doc me
    const CouplingContext& couplingContext() const
    {
        if (!couplingContext_)
            DUNE_THROW(Dune::InvalidStateException, "Coupling context not available!");
        return *couplingContext_;
    }

private:
    const ElementSolution* elemSol_ = nullptr;
    const ElementVariables* elemVars_ = nullptr;
    const CouplingContext* couplingContext_ = nullptr;
    std::optional<TimeLevel> timeLevel_;
};

} // end namespace Dumux

#endif
