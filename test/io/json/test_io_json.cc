// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test json
 */

#include <iostream>

#include <dumux/io/json.hh>

#include <dune/common/fvector.hh>

int main()
{
    using namespace Dumux::Json;

    // check conversion to FieldVector
    {
        JsonTree j = R"({
            "array": [1, 2, 3]
        })"_json;

        // convert to Dune::FieldVector
        auto v = j["array"].template get<Dune::FieldVector<double, 3>>();

        // convert back to json
        j["array"] = v;
    }

    std::cout << "All json tests passed." << std::endl;

    return 0;
}
