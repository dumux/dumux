// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Particles
 * \brief A simple particle
 */
#ifndef DUMUX_PARTICLES_PARTICLE_HH
#define DUMUX_PARTICLES_PARTICLE_HH

#include <dune/common/fvector.hh>

namespace Dumux::Particles {

/*!
 * \ingroup Particles
 * \brief a basic particle
 */
template<int dimWorld, class ctype = double>
class Particle
{
public:
    static constexpr int coordDimension = dimWorld;
    using CoordScalar = ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    Particle(std::size_t id, const GlobalPosition& p, bool active = true)
    : id_(id), position_(p), active_(active) {}

    Particle(std::size_t id, bool active = true)
    : Particle(id, GlobalPosition(0.0), active) {}

    const GlobalPosition& position() const
    { return position_; }

    std::size_t id() const
    { return id_; }

    bool isActive() const
    { return active_; }

    //! set a new global position
    void setPosition(const GlobalPosition& pos)
    { position_ = pos; }

    //! move by length in direction
    void move(const CoordScalar length, const GlobalPosition& direction)
    { position_.axpy(length, direction); }

    void deactivate()
    { active_ = false; }

    void activate()
    { active_ = true; }

private:
    std::size_t id_;
    GlobalPosition position_;
    bool active_;
};

} // end namespace Dumux::Particles

#endif
