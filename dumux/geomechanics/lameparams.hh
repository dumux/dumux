// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsModels
 * \brief \copydoc Dumux::LameParams
 */
#ifndef DUMUX_GEOMECHANICS_LAME_PARAMS_HH
#define DUMUX_GEOMECHANICS_LAME_PARAMS_HH

namespace Dumux {

/*!
 * \ingroup GeomechanicsModels
 * \brief Structure encapsulating the lame parameters
 */
template<class Scalar>
struct LameParams
{
    //! Default constructor
    LameParams() = default;

    //! Constructor taking lambda and mu directly
    LameParams(Scalar lambda, Scalar mu)
    : lambda_(lambda) , mu_(mu)
    {}

    //! Return the first lame parameter
    Scalar lambda() const
    { return lambda_; }

    //! Return the second lame parameter
    Scalar mu() const
    { return mu_; }

    //! set the first lame parameter
    void setLambda(Scalar lambda)
    { lambda_ = lambda; }

    //! set the second lame parameter
    void setMu(Scalar mu)
    { mu_ = mu; }

private:
    Scalar lambda_;
    Scalar mu_;
};
} // end namespace Dumux
#endif
