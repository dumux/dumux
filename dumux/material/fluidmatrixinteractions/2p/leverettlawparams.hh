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
 * \ingroup Fluidmatrixinteractions
 * \brief Specification of the material parameters
 *       for the van Genuchten constitutive relations.
 */
#ifndef LEVERETT_LAW_PARAMS_HH
#define LEVERETT_LAW_PARAMS_HH

#include <dune/common/float_cmp.hh>

namespace Dumux
{
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Specification of the material parameters
 *       for the van Genuchten constitutive relations.
 *
 *       In this implementation setting either the \f$\mathrm{n}\f$ or \f$\mathrm{m}\f$ shape parameter
 *       automatically calculates the other. I.e. they cannot be set independently.
 */
template<class BaseParams>
class LeverettLawParams : public BaseParams
{
public:
    using Scalar = typename BaseParams::Scalar;

    LeverettLawParams()
    {
                referencePorosity_     = getParam<Scalar>("SpatialParams.referencePorosity", 0.11);
        referencePermeability_ = getParam<Scalar>("SpatialParams.referencePermeability", 2.23e-14);
    }

    LeverettLawParams(Scalar l)
    {
        setLeverettFactor(l);
    }


    void setLeverettFactor(Scalar l)
    {leverettFactor_ = l;}

    void setLeverettFactor(Scalar porosity, Scalar permeability)
    {
        leverettFactor_ = pow(referencePermeability_/permeability*porosity/referencePorosity_, 0.5);
    }


    Scalar leverettFactor() const
    {return leverettFactor_; }

private:
    Scalar referencePorosity_;
    Scalar referencePermeability_;
    Scalar leverettFactor_;
};
} // namespace Dumux

#endif
