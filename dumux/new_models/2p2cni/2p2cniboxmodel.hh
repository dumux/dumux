/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_NEW_2P2CNI_BOX_MODEL_HH
#define DUMUX_NEW_2P2CNI_BOX_MODEL_HH

#include <dumux/new_models/2p2c/2p2cboxmodel.hh>

#include <dumux/new_models/2p2cni/2p2cniboxjacobian.hh>
#include <dumux/new_models/2p2cni/2p2cniproperties.hh>

namespace Dune
{
/**
 * \brief Non-isothermal two-phase two-component model with Pw and
 *        Sn/X as primary unknowns.
 *
 * This implements a non-isothermal two-phase two-component model with
 * Pw and Sn/X as primary unknowns. You can use Pn and Sw/X as primary
 * variables if you set the Formulation property to pNsW.
 */
template<class TypeTag>
class TwoPTwoCNIBoxModel
    : public TwoPTwoCBoxModelBase<TypeTag, TwoPTwoCNIBoxModel<TypeTag> >
{
    typedef TwoPTwoCNIBoxModel<TypeTag>                           ThisType;
    typedef TwoPTwoCBoxModelBase<TypeTag, ThisType>               ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))       Problem;

public:
    TwoPTwoCNIBoxModel(Problem &prob)
        : ParentType(prob)
    {
    }
};

}

#endif
