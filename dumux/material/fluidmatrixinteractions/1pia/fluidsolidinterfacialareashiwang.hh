// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief Description of a interfacial area between solid and fluid phase
 */
#ifndef FLUIDSOLID_INTERFACIALAREA_SHI_WANG_HH
#define FLUIDSOLID_INTERFACIALAREA_SHI_WANG_HH

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Description of a interfacial area between solid and fluid phase
 */
template<class Scalar>
class FluidSolidInterfacialAreaShiWang
{
public:
    /*!
     * \brief Relation for the interfacial area between a fluid and a solid phase
     * after Shi & Wang, Transport in porous media (2011)
     *
     * \return interfacial area
     */
    static Scalar fluidSolidInterfacialArea(const Scalar porosity,
                                            const Scalar characteristicLength)
    { return 6.0*(1.0-porosity)/characteristicLength; }
};

} // end namespace Dumux

#endif
