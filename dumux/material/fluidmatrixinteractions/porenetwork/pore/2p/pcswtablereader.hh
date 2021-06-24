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
 * \ingroup Components
 * \brief A generic template for tabulated material laws that depend
 *        on one parameter.
 */
#ifndef DUMUX_TABULATED_PCSW_HH
#define DUMUX_TABULATED_PCSW_HH

#include <dune/common/float_cmp.hh>
#include <dumux/common/exceptions.hh>

namespace Dumux {
/*!
 * \ingroup Components
 * \brief A generic template for tabulated material laws that depend
 *        on one parameter.
 */
template <class Traits>
class TabulatedPcSwProperties
{
    using Scalar = typename Traits::Scalar;
    enum { numSteps = Traits::numSteps};

public:
    // user default constructor (we can't use "= default" here to satisfy older clang compilers since this class is used as a static data member)
    TabulatedPcSwProperties() {}
// sw functions
    Scalar minSw() const
    { return Traits::minSw; }

    Scalar maxSw() const
    { return Traits::maxSw; }

    bool applies(Scalar sw) const
    {
        return minSw() <= sw && sw <= maxSw();
    }

    Scalar at(Scalar pc) const //get sw corresponding to pc
    {
        if (!appliesPc(pc))
        {
            if (pc<minPc())
                pc=minPc();
            if(pc>maxPc())
                pc=maxPc();
        }

        int i = findPcIdx_(pc);
        return val(i);
        /*Scalar swAtI = swAt_(i);
        Scalar swAtI1 = swAt_(i + 1);

        Scalar alpha = (sw - swAtI)/(swAtI1 - swAtI);
        return alpha;*/
    }

    Scalar val(int i) const
    {
#if !defined NDEBUG
        if (i < 0 || i >= Traits::numSteps) {
            DUNE_THROW(NumericalProblem,
                       "Attempt to access element ("
                       << i << ") on a " << Traits::name << " table of size ("
                       << Traits::numSteps << ")\n");
        }
#endif
        return Traits::vals[i];
    }

protected:
    int findSwIdx_(Scalar sw) const
    {
        if (Dune::FloatCmp::eq<Scalar>(sw, maxSw()))
            return numSteps ;
        //const int result = static_cast<int>((sw - minSw())/(maxSw() - minSw())*(numSteps - 1));
        int result = 0.0;
        if (sw > 0.9995905)
            result = static_cast<int>((sw - 0.9995905)/(maxSw() - 0.9995905)*(2444 - 1));
        else
            result = static_cast<int>((sw - minSw())/(0.9995905-minSw())*(614 - 1));
        using std::min;
        using std::max;
        return max(0, min(result, numSteps - 2));
    }

    Scalar swAt_(int i) const
    { return i*(maxSw() - minSw())/(numSteps - 1) + minSw(); }

public:
// pc functions
    Scalar minPc() const
    { return Traits::minPc; }

    Scalar maxPc() const
    { return Traits::maxPc; }

    bool appliesPc(Scalar pc) const
    {
        return minPc() <= pc && pc <= maxPc();
    }

    Scalar atSw(Scalar sw) const //get pc corresponding to sw
    {
        if (!applies(sw))
        {
            if (sw<minSw())
                sw=minSw();
            if(sw>maxSw())
                sw=maxSw();
        }

        int i = findSwIdx_(sw);
        return val(i);
        /*Scalar pcAtI = pcAt_(i);
        Scalar pcAtI1 = pcAt_(i + 1);

        Scalar alpha = (pc - pcAtI)/(pcAtI1 - pcAtI);
        return alpha;*/
    }

    Scalar valPc(int i) const
    {
#if !defined NDEBUG
        if (i < 0 || i >= Traits::numSteps) {
            DUNE_THROW(NumericalProblem,
                       "Attempt to access element ("
                       << i << ") on a " << Traits::name << " table of size ("
                       << Traits::numSteps << ")\n");
        }
#endif
        return Traits::vals[i];
    }

protected:
    int findPcIdx_(Scalar pc) const
    {
        if (Dune::FloatCmp::eq<Scalar>(pc, maxPc()))
            return numSteps - 2;
        const int result = static_cast<int>((pc - minPc())/(maxPc() - minPc())*(numSteps - 1));

        using std::min;
        using std::max;
        return max(0, min(result, numSteps - 2));
    }

    //Scalar pcAt_(int i) const
    //{ return i*(maxPc() - minPc())/(numSteps - 1) + minPc(); }
};
} // end namespace Dumux

#endif
