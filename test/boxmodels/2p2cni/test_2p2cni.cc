// $Id: test_2p2cni.cc 3779 2010-06-24 07:07:56Z bernd $
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Melanie Darcis                               *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
#include "config.h"
#include "waterairproblem.hh"
#include <dumux/common/start.hh>

int main(int argc, char** argv)
{
    typedef TTAG(WaterAirProblem) ProblemTypeTag;
    return Dumux::startFromDGF<ProblemTypeTag>(argc, argv);
}
