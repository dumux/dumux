// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!\file
 * \ingroup Common
 * \brief Adapter to expose a multi-DOF mapper interface for single- and multi-DOF mappers.
 */
#ifndef DUMUX_COMMON_MULTIMAPPERVIEW_HH
#define DUMUX_COMMON_MULTIMAPPERVIEW_HH

#include <array>

namespace Dumux::Detail {

template<class Mapper>
class MultiMapperView
{
public:
    explicit constexpr MultiMapperView(const Mapper& mapper)
    : mapper_(mapper)
    {}

    template<class Entity>
    decltype(auto) index(const Entity& entity) const
    { return mapper_.index(entity); }

    template<class Entity>
    auto indices(const Entity& entity) const
    {
        if constexpr (requires { mapper_.indices(entity); })
            return mapper_.indices(entity);
        else
            return std::array<decltype(mapper_.index(entity)), 1>{ mapper_.index(entity) };
    }

private:
    const Mapper& mapper_;
};

} // end namespace Dumux::Detail

namespace Dumux {

template<class Mapper>
constexpr auto asMultiMapper(const Mapper& mapper)
{
    return Detail::MultiMapperView<Mapper>(mapper);
}

} // end namespace Dumux

#endif
