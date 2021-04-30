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
 * \ingroup Common
 * \brief Test for the SolutionVector wrapper class
 */
#include <config.h>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/common/partial.hh>
#include <dumux/common/solutionvector.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>


int main(int argc, char** argv)
{
    using namespace Dumux;

    using StandardPrivars = Dune::FieldVector<double, 2>;
    using SwitchablePrimaryVariables = Dumux::SwitchablePrimaryVariables<StandardPrivars, int>;

    using StandardSolutionVector = Dumux::SolutionVector<StandardPrivars>;
    using SolutionVectorWithSwitch = Dumux::SolutionVector<SwitchablePrimaryVariables>;
    using MultiTypeSolutionVector = Dumux::SolutionVector<StandardPrivars, StandardPrivars>;
    using MultiTypeSolutionVectorWithSwitch = Dumux::SolutionVector<StandardPrivars, SwitchablePrimaryVariables>;


    StandardSolutionVector xStandard(10);
    xStandard = 123.0;
    static_assert(std::is_same_v<decltype(xStandard), Dune::BlockVector<StandardPrivars>>);

    static_assert(std::is_same_v<Dumux::Istl::BlockVectorWithState<Dune::BlockVector<StandardPrivars>,
                                                                   Dumux::SwitchablePrimaryVariables<StandardPrivars, int>>,
                  SolutionVectorWithSwitch>);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////// SolutionVectorWithSwitch
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    SolutionVectorWithSwitch xWithSwitch(10);
    xWithSwitch = 1.0;

    // assigen a "real" PriVars object to the view
    SwitchablePrimaryVariables pv(123.0);
    pv.setState(1);
    xWithSwitch[1] = pv;
    if (xWithSwitch[1][0] != 123.0)
        DUNE_THROW(Dune::InvalidStateException, "Failed");
    if (xWithSwitch[1].state() != 1)
        DUNE_THROW(Dune::InvalidStateException, "Failed");

    // assign the view to a "real" PriVars object
    xWithSwitch[0] = 222.0;
    xWithSwitch[0].setState(2);
    pv = xWithSwitch[0];
    if (pv[0] != 222.0)
        DUNE_THROW(Dune::InvalidStateException, "Failed");
    if (pv.state() != 2)
        DUNE_THROW(Dune::InvalidStateException, "Failed");

    // test reference
    auto& xWithSwitchRef = xWithSwitch;
    xWithSwitch[1] *= 2.0;
    if (xWithSwitchRef[1][0] != 246.0)
        DUNE_THROW(Dune::InvalidStateException, "Failed");

    // test const object
    const auto xWithSwitchConst = xWithSwitch;
    if (xWithSwitchConst[0][0] != 222.0)
        DUNE_THROW(Dune::InvalidStateException, "Failed");

    // test native type
    auto xWithSwitchNative = native(xWithSwitch);
    static_assert(std::is_same_v<decltype(xWithSwitchNative), Dune::BlockVector<Dune::FieldVector<double, 2>>>);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////// MultiTypeSolutionVectorWithSwitch
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    MultiTypeSolutionVectorWithSwitch xMultiWithSwitch;

    // assing other solution vector
    xMultiWithSwitch[Dune::index_constant<0>{}] = xStandard;
    xMultiWithSwitch[Dune::index_constant<1>{}] = xWithSwitch;

    // iterate over vector
    for (auto& x : xMultiWithSwitch[Dune::index_constant<1>{}])
    {
        x = 333.0;
        static_assert(std::is_same_v<std::decay_t<decltype(x)>, Dune::FieldVector<double, 2>>);
        // TODO: maybe change this such that x is a View object, too (Dumux::BlockVector needs to be changed)
    }

    if (xMultiWithSwitch[Dune::index_constant<1>{}][0][0] != 333.0)
        DUNE_THROW(Dune::InvalidStateException, "Failed");

    // changing the view changes the underlying MultiTypeSolutionVectorWithSwitch, too
    // TODO: is this what we want?
    auto view = xMultiWithSwitch[Dune::index_constant<1>{}];
    view = 555.0;
    if (xMultiWithSwitch[Dune::index_constant<1>{}][0][0] != 555.0)
         DUNE_THROW(Dune::InvalidStateException, "Failed");

    // check native type
    static_assert(std::is_same_v<std::decay_t<decltype(native(xMultiWithSwitch))>,
                                 Dune::MultiTypeBlockVector<Dune::BlockVector<Dune::FieldVector<double, 2>>,
                                                            Dune::BlockVector<Dune::FieldVector<double, 2>>>>);

    // test "simpler" MultiTypeSolutionVector without SwitchablePrimaryVariables // TODO [] operator will behave differenty, not yielding a view (yet)
    // but just forwarding to une::MultiTypeBlockVector's []
    static_assert(std::is_same_v<typename MultiTypeSolutionVector::NativeType,
                                 Dune::MultiTypeBlockVector<Dune::BlockVector<Dune::FieldVector<double, 2>>,
                                                            Dune::BlockVector<Dune::FieldVector<double, 2>>>>);

    // TODO implement/check partial

} // end main
