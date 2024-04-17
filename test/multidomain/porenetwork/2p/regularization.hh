// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief This file contains functions for regularizing throat conductivities
 */
#ifndef DUMUX_PNM_THROAT_REGULARIZATION_HH
#define DUMUX_PNM_THROAT_REGULARIZATION_HH

#include <string>
#include <cmath>
#include <numeric>

#include <dune/common/exceptions.hh>
#include <dumux/common/spline.hh>

namespace Dumux::PoreNetwork::Throat {

//! Collection of different regularization functions
enum class RegFunction
{ sine, linear, spline };

//! Get the function name from a string description
inline std::string functionToString(RegFunction s)
{
    switch (s)
    {
        case RegFunction::sine: return "Sine";
        case RegFunction::linear: return "Linear";
        case RegFunction::spline: return "Spline";
        default: DUNE_THROW(Dune::InvalidStateException, "Unknown function!");
    }
}

//! Get the function from its name (string)
inline RegFunction functionFromString(const std::string& s)
{
    if (s == functionToString(RegFunction::sine)) return RegFunction::sine;
    else if (s == functionToString(RegFunction::linear)) return RegFunction::linear;
    else if (s == functionToString(RegFunction::spline)) return RegFunction::spline;
    else DUNE_THROW(Dune::InvalidStateException, s << " is not a valid function name");
}

template <class Scalar>
class Regularization
{

public:
    Regularization(std::string s, Scalar delta) : delta_(delta)
    {
        mode_ = functionFromString(s);
    }

    Scalar eval(Scalar v, bool invaded) const
    {
        using std::min; using std::max;
        if (mode_ == RegFunction::sine)
        {
            // TODO Also only use delta and not 2delta interval
            using std::sin;
            if(!invaded)
            {
                using std::sin; using std::min; using std::max;
                v = max(0.0,min(v,delta_));
                return 0.5*(1+sin(M_PI*(v-0.5*delta_)/(delta_)));
            }
            else
            {
                using std::sin; using std::min; using std::max;
                v = max(-delta_,min(v,0.0));
                return 0.5*(1+sin(M_PI*(v+0.5*delta_)/(delta_)));
            }
        }
        else if(mode_ == RegFunction::linear)
        {
            if(!invaded)
            {
                v = max(0.0,min(v,delta_));
                return v/delta_;
            }
            else
            {
                v = max(-delta_,min(v,0.0));
                return 1.0 + v/delta_;
            }
        }
        else if(mode_ == RegFunction::spline)
        {
            // TODO do not assume zero derivatives but generalize this
            if(!invaded)
            {
                auto spline = Spline<Scalar>(0, delta_, // x0, x1
                                             0.0, 1.0, // y0, y1
                                             0.0, 0.0); // deriv0, deriv1
                v = max(0.0,min(v,delta_));
                return spline.eval(v);
            }
            else
            {
                auto spline = Spline<Scalar>(-delta_, 0, // x0, x1
                                             0.0, 1.0, // y0, y1
                                             0.0, 0.0); // deriv0, deriv1
                v = max(-delta_,min(v,0.0));
                return spline.eval(v);
            }
        }
        else DUNE_THROW(Dune::InvalidStateException, "Unkown regularization function.");
    }

private:
    Scalar delta_;
    RegFunction mode_;
};


} // end Dumux::PoreNetwork::ThroatRegularization

#endif
