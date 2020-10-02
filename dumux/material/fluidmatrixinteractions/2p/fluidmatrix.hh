// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \brief   Implementation helper for capillary pressure and
 *          relative permeability <-> saturation relations for two-phase models
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_FLUIDMATRIX_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_FLUIDMATRIX_HH



namespace Dumux::FluidMatrix {

template<class T, class StoreT = const T&>
struct PcSw
{
    using BasicParams = typename T::BasicParams;
    using EffToAbsParams = typename T::EffToAbsParams;

    PcSw(const T& i) : impl_(i) {}
    const T& pcSw() const { return impl_; }
private:
    StoreT impl_;
};

template<class T, class StoreT = const T&>
struct KrSw
{
    using BasicParams = typename T::BasicParams;
    using EffToAbsParams = typename T::EffToAbsParams;

    KrSw(const T& i) : impl_(i) {}
    const T& krSw() const { return impl_; }
private:
    StoreT impl_;
};

// template<class T, class StoreT = const T&>
// struct KrSw
// {
//     using PcSw = T;
//     PcSw(const T& i) : impl_(i) {}
//     const T& pcSw() const { return impl_; }
// private:
//     StoreT impl_;
// };

// template<class T, class StoreT = const T&>
// struct InterfacialAreaSwPc
// {
//     using PcSw = T;
//     PcSw(const T& i) : impl_(i) {}
//     const T& pcSw() const { return impl_; }
// private:
//     StoreT impl_;
// };


template<class... Laws>
struct FluidMatrixInteraction : public Laws...
{
    template<class... Args>
    FluidMatrixInteraction(Args&&... args) : Laws(std::forward<Args>(args))... {}
};

} // end namespace Dumux::FluidMatrix

#endif
