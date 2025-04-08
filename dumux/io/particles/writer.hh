// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Particles
 * \brief A particle output writer
 */

#ifndef DUMUX_IO_PARTICLES_WRITER_HH
#define DUMUX_IO_PARTICLES_WRITER_HH

#include <memory>
#include <fstream>
#include <iomanip>
#include <vector>
#include <ranges>

#include <dumux/io/gridwriter.hh>

#ifdef DUMUX_HAVE_GRIDFORMAT
#include <gridformat/gridformat.hpp>

namespace Dumux::IO::Particles {

/**
 * \file
 * \ingroup Particles
 * \brief A time series writer for particles
 */
template <class ParticleCloud>
class Writer
{
public:
    template <class Format>
    Writer(Format fmt, std::shared_ptr<ParticleCloud> cloud, const std::string& name)
    : cloud_(std::move(cloud))
    , name_(name)
    , writer_(std::move(fmt), *cloud_, name)
    {}

    Writer(std::shared_ptr<ParticleCloud> cloud, const std::string& name)
    : Writer(
        IO::Format::pvd_with(IO::Format::vtp.with({.encoder = IO::Encoding::raw})
    ), std::move(cloud), name)
    {}

    template<std::floating_point T>
    std::string write(T time)
    { return writer_.write(time); }

    template<Detail::Container C, GridFormat::Concepts::Scalar T>
    void addParticleData(
        const C& data, const std::string& name,
        const GridFormat::Precision<T>& prec
    ){
        writer_.set_cell_field(name, [&] (const auto& p) -> std::ranges::range_value_t<C> {
            return data[p.id()];
        }, prec);
    }

    template<Detail::Container C>
    void addParticleData(
        const C& data, const std::string& name
    ){
        this->addParticleData(
            data, name,
            GridFormat::Precision<GridFormat::MDRangeScalar<C>>{}
        );
    }

private:
    std::shared_ptr<ParticleCloud> cloud_;
    const std::string name_;
    GridFormat::Writer<std::remove_cvref_t<ParticleCloud>> writer_;
};

} // end namespace Dumux::IO::Particles

#endif // DUMUX_HAVE_GRIDFORMAT

#endif
