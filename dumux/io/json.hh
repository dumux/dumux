// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Collection of json classes from JSON for Modern C++ library
 */

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dumux/io/json/json.hpp>

namespace Dumux::Json {
using JsonTree = nlohmann::json;
} // end namespace Dumux::Json

namespace nlohmann {

template <typename ctype, int k>
struct adl_serializer<Dune::FieldVector<ctype, k>>
{
    static void to_json(json& j, const Dune::FieldVector<ctype, k>& fv)
    {
        j = json::array();
        for (int i = 0; i < k; ++i) {
            j.push_back(fv[i]);
        }
    }

    static void from_json(const json& j, Dune::FieldVector<ctype, k>& fv)
    {
        if (!j.is_array())
            DUNE_THROW(Dune::IOError, "json: Cannot convert to FieldVector, not an array");

        if (j.size() != k)
            DUNE_THROW(Dune::IOError, "json: Cannot convert to FieldVector of size " << k << " from array of size " << j.size());

        for (int i = 0; i < k; ++i)
            fv[i] = j[i].template get<ctype>();
    }
};

} // end namespace nlohmann
