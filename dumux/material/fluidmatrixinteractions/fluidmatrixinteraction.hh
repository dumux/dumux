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
 * \brief   Wrapper type to combine an arbitrary number of different laws
 *          for fluid-matrix interaction (e.g., pc-Sw-curves).
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_INTERACTIONS_FLUIDMATRIX_INTERACTION_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_INTERACTIONS_FLUIDMATRIX_INTERACTION_HH

#include <type_traits>
#include <utility>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Wrapper type to combine an arbitrary number of different laws
 *        for fluid-matrix interaction (e.g., pc-Sw-curves).
 */
template<class... Laws>
struct FluidMatrixInteraction : public Laws...
{
   FluidMatrixInteraction(Laws&&... laws) : Laws(std::forward<Laws>(laws))... {}
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Helper function to create an FluidMatrixInteraction object containing an arbitrary number
 *        of fluid matrix interaction laws (e.g., pc-Sw curves and interfacial area laws).
 *        To be used in the spatial parameters.
 */
template<class... Laws>
auto makeFluidMatrixInteraction(Laws&&... laws)
{
    return FluidMatrixInteraction(wrap(std::forward<Laws>(laws))...);
}

} // end namespace Dumux

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Adapter to inherit from, allowing the inheriting class to be wrapped
 *        by the @ref makeFluidMatrixInteraction function.
 */
template<class A, template<class> class Wrapper>
struct Adapter
{
    template<class T, std::enable_if_t<std::is_same_v<A, std::decay_t<T>>, int> = 0>
    friend auto wrap(T&& t) { return Wrapper(std::forward<T>(t)); }
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Wrapper type for laws providing pc-Sw and kr-Sw rules.
 */
template<class T>
class PcKrSw
{
public:
    using Scalar = typename std::decay_t<T>::Scalar;

    using PcKrSwType = T;

    PcKrSw(T&& impl) : impl_(std::forward<T>(impl)) {}

    Scalar pc(const Scalar sw) const { return impl_.pc(sw); }
    Scalar dpc_dsw(const Scalar sw) const { return impl_.dpc_dsw(sw); }
    Scalar endPointPc() const { return impl_.endPointPc(); }
    Scalar sw(const Scalar pc) const { return impl_.sw(pc); }
    Scalar dsw_dpc(const Scalar pc) const { return impl_.dsw_dpc(pc); }
    Scalar krw(const Scalar sw) const { return impl_.krw(sw); }
    Scalar dkrw_dsw(const Scalar sw) const { return impl_.dkrw_dsw(sw); }
    Scalar krn(const Scalar sw) const { return impl_.krn(sw); }
    Scalar dkrn_dsw(const Scalar sw) const { return impl_.dkrn_dsw(sw); }

    const T& pcSwCurve() const { return impl_; }
    const T& krSwCurve()  const { return impl_; }

private:
    T impl_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Deduction guide for the PcKrSw class.
 *        Makes sure that PcKrSw stores a copy of T if
 *        the constructor is called with a temporary object.
 */
template<typename T>
PcKrSw(T&&) -> PcKrSw<T>;

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Wrapper type for multiphase interface laws providing pc-S and kr-S rules.
 */
template<class T>
class MultiPhasePcKrSw
{
public:
    using Scalar = typename std::decay_t<T>::Scalar;

    MultiPhasePcKrSw(T&& impl) : impl_(std::forward<T>(impl)) {}

    template<class FS>
    auto capillaryPressures(const FS& fs, int wp) const { return impl_.capillaryPressures(fs, wp); }
    template<class FS>
    auto relativePermeabilities(const FS& fs, int wp) const { return impl_.relativePermeabilities(fs, wp); }

    const T& multiPhasePcKrS() const { return impl_; }

private:
    T impl_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Deduction guide for the MultiPhasePcKrSw class.
 *        Makes sure that MultiPhasePcKrSw stores a copy of T if
 *        the constructor is called with a temporary object.
 */
template<typename T>
MultiPhasePcKrSw(T&&) -> MultiPhasePcKrSw<T>;

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Wrapper type for 3p interface laws providing pc-S and kr-S rules.
 */
template<class T>
struct ThreePhasePcKrSw
{
    using Scalar = typename std::decay_t<T>::Scalar;
    using value_type = T;

    using PcKrSwType = T;

    ThreePhasePcKrSw(T&& impl) : impl_(std::forward<T>(impl)) {}

    Scalar pcgw(const Scalar sw, const Scalar sn) const { return impl_.pcgw(sw, sn); }
    Scalar pcnw(const Scalar sw, const Scalar sn) const { return impl_.pcnw(sw, sn); }
    Scalar pcgn(const Scalar sw, const Scalar sn) const { return impl_.pcgn(sw, sn); }
    Scalar pcAlpha(const Scalar sw, const Scalar sn) const { return impl_.pcAlpha(sw, sn); }

    Scalar krw(const Scalar sw, const Scalar sn) const { return impl_.krw(sw, sn); }
    Scalar krn(const Scalar sw, const Scalar sn) const { return impl_.krn(sw, sn); }
    Scalar krg(const Scalar sw, const Scalar sn) const { return impl_.krn(sw, sn); }
    Scalar kr(const int phaseIdx, const Scalar sw, const Scalar sn) const { return impl_.kr(phaseIdx, sw, sn); }

    const T& pcSwCurve() const { return impl_; }
    const T& krSwCurve()  const { return impl_; }
private:
    T impl_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Deduction guide for the ThreePhasePcKrSw class.
 *        Makes sure that ThreePhasePcKrSw stores a copy of T if
 *        the constructor is called with a temporary object.
 */
template<typename T>
ThreePhasePcKrSw(T&&) -> ThreePhasePcKrSw<T>;

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Wrapper type for laws providing rules for the wetting-nonwetting interfacial area.
 */
template<class T>
class WettingNonwettingInterfacialAreaPcSw
{
public:
    WettingNonwettingInterfacialAreaPcSw(T&& impl) : impl_(std::forward<T>(impl)) {}
    const T& wettingNonwettingInterface() const { return impl_; }
private:
    T impl_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Deduction guide for the WettingNonwettingInterfacialAreaPcSw class.
 *        Makes sure that WettingNonwettingInterfacialAreaPcSw stores a copy of T if
 *        the constructor is called with a temporary object.
 */
template<typename T>
WettingNonwettingInterfacialAreaPcSw(T&&) -> WettingNonwettingInterfacialAreaPcSw<T>;

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Wrapper type for laws providing rules for the wetting-solid interfacial area.
 */
template<class T>
class WettingSolidInterfacialAreaPcSw
{
public:
    WettingSolidInterfacialAreaPcSw(T&& impl) : impl_(std::forward<T>(impl)) {}
    const T& wettingSolidInterface() const { return impl_; }
private:
    T impl_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Deduction guide for the WettingSolidInterfacialAreaPcSw class.
 *        Makes sure that WettingSolidInterfacialAreaPcSw stores a copy of T if
 *        the constructor is called with a temporary object.
 */
template<typename T>
WettingSolidInterfacialAreaPcSw(T&&) -> WettingSolidInterfacialAreaPcSw<T>;

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Wrapper type for laws providing rules for the nonwetting-solid interfacial area.
 */
template<class T>
class NonwettingSolidInterfacialAreaPcSw
{
public:
    NonwettingSolidInterfacialAreaPcSw(T&& impl) : impl_(std::forward<T>(impl)) {}
    const T& nonwettingSolidInterface() const { return impl_; }
private:
    T impl_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Deduction guide for the NonwettingSolidInterfacialAreaPcSw class.
 *        Makes sure that NonwettingSolidInterfacialAreaPcSw stores a copy of T if
 *        the constructor is called with a temporary object.
 */
template<typename T>
NonwettingSolidInterfacialAreaPcSw(T&&) -> NonwettingSolidInterfacialAreaPcSw<T>;


/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Wrapper type for adsorption laws.
 */
template<class T>
class Adsorption
{
public:
    using value_type = T;

    Adsorption(T&& impl) : impl_(std::forward<T>(impl)) {}
    const T& adsorptionModel() const { return impl_; }
private:
    T impl_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Deduction guide for the Adsorption class.
 *        Makes sure that Adsorption stores a copy of T if
 *        the constructor is called with a temporary object.
 */
template<typename T>
Adsorption(T&&) -> Adsorption<T>;

}  // end namespace Dumux::FluidMatrix

#endif
