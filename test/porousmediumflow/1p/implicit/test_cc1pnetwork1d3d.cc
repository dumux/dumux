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
 *
 * \brief test for the one-phase CC model
 */
#include <config.h>
#include "tubesproblem.hh"
#include <dumux/common/start.hh>

int main(int argc, char** argv)
{
#if HAVE_DUNE_FOAMGRID
    typedef TTAG(TubesTestCCTpfaProblem) ProblemTypeTag;
    return Dumux::start<ProblemTypeTag>(argc, argv, [](const char *, const std::string &){});
#else
#warning External grid module dune-foamgrid needed to run this example.
    std::cerr << "Test skipped, it needs dune-foamgrid!" << std::endl;
    return 77;
#endif
}
