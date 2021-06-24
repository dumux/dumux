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
 * \brief Base classes for standard pore-local pc-Sw curves.
 */
#ifndef DUMUX_PNM_2P_DROP_LOCAL_RULES_HH
#define DUMUX_PNM_2P_DROP_LOCAL_RULES_HH
#include <cmath>
#include <dumux/common/parameters.hh>
#include <dumux/common/optionalscalar.hh>
#include <dumux/porenetwork/common/poreproperties.hh>
//#include <dumux/material/fluidmatrixinteractions/2p/noregularization.hh>
//#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/common/typetraits/typetraits.hh>
#include "singleshapelocalrules.hh"
#include "pcswdroptables.hh"
#include "pcswtablereader.hh"


namespace Dumux::PoreNetwork::FluidMatrix {

/*!
 * \brief The parameter type
 * \tparam Scalar The scalar type
 */
template<class Scalar>
struct DropParams
{
    DropParams() = default;

    DropParams& setPoreInscribedRadius(Scalar r) { radius_ = r; return *this;}
    DropParams& setPoreShape(Pore::Shape s) { shape_ = s; return *this;}
    DropParams& setSurfaceTension(Scalar st) { surfaceTension_ = st; return *this;}
    DropParams& setThroatInscribedRadius(Scalar r) { radiusThroat_ = r; return *this;}
    DropParams& setCriticalContactAngle(Scalar theta) {thetaCrit_ = theta; return *this;}

    template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
    DropParams(const SpatialParams& spatialParams,
                       const Element& element,
                       const SubControlVolume& scv,
                       const ElemSol& elemSol)
    : shape_(spatialParams.gridGeometry().poreGeometry()[scv.dofIndex()])
    , radius_(spatialParams.poreInscribedRadius(element, scv, elemSol))
    , radiusThroat_(spatialParams.throatInscribedRadius(element))
    {
        static const Scalar surfaceTension = getParam<Scalar>("SpatialParams.SurfaceTension", 0.0725); // TODO
        surfaceTension_ = surfaceTension;
        static const Scalar thetaCrit = getParam<Scalar>("SpatialParams.CriticalContactAngle", (90-120)/360*2*M_PI); // TODO
        thetaCrit_ = thetaCrit;
    }

    template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
    void update(const SpatialParams& spatialParams,
                const Element& element,
                const SubControlVolume& scv,
                const ElemSol& elemSol)
    {
        const auto& gridGeometry = spatialParams.gridGeometry();
        shape_ = gridGeometry.poreGeometry()[scv.dofIndex()];
        radius_ = spatialParams.poreInscribedRadius(element, scv, elemSol);
        radiusThroat_ = spatialParams.throatInscribedRadius(element);

        static const Scalar surfaceTension = getParam<Scalar>("SpatialParams.SurfaceTension", 0.0725); // TODO
        surfaceTension_ = surfaceTension;

        static const Scalar thetaCrit = getParam<Scalar>("SpatialParams.CriticalContactAngle", (90-120)/360*2*M_PI); // TODO
        thetaCrit_ = thetaCrit;
    }

    Pore::Shape poreShape() const { return shape_; }

    Scalar poreInscribedRadius() const { return radius_; }

    Scalar surfaceTension() const { return surfaceTension_; }
    Scalar throatInscribedRadius() const { return radiusThroat_; }
    Scalar criticalContactAngle() const { return thetaCrit_; }

    bool operator== (const DropParams& p) const
    {
        return Dune::FloatCmp::eq(radius_, p.radius_, 1e-6)
               && Dune::FloatCmp::eq(surfaceTension_, p.surfaceTension_, 1e-6)
               && shape_ == p.shape_
               && Dune::FloatCmp::eq(radiusThroat_, p.radiusThroat_, 1e-6)
               && Dune::FloatCmp::eq(thetaCrit_, p.thetaCrit_, 1e-6);
    }

private:
    Pore::Shape shape_;
    Scalar radius_;
    Scalar surfaceTension_;
    Scalar radiusThroat_;
    Scalar thetaCrit_;
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of tabulated capillary pressure-saturation curve for drop pores
 */
template<Pore::Shape shape>
struct DropLocalRules
{
    template<class Scalar>
    using Params = DropParams<Scalar>;

    template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
    static auto makeParams(const SpatialParams& spatialParams,
                           const Element& element,
                           const SubControlVolume& scv,
                           const ElemSol& elemSol)
    {
        using Scalar = std::decay_t<decltype(spatialParams.poreInscribedRadius(element, scv, elemSol))>;
        return Params<Scalar>(spatialParams, element, scv, elemSol);
    }

/*!
     * \brief The capillary pressure-saturation curve
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{S_{w,i}}\f$ at pore \f$i\f$
     * \param params The parameters container
     */
    template<class Scalar>
    //template<class PcSwDropTables>
    static Scalar pc(Scalar sw, const Params<Scalar>& params)
    {
        if (sw < swCrit(params))
        {
            const Scalar sigma = params.surfaceTension();
            const Scalar throatRadius = params.throatInscribedRadius();
            const Scalar throatArea = M_PI*throatRadius*throatRadius;
            const Scalar theta = asin(std::sqrt(throatArea/(M_PI*dropRadius)));

            return 2*sigma*std::cos(theta)/throatRadius;
        }
        PcSwTables::TabulatedPc tabulatedPc;
        return tabulatedPc.TabulatedPcSwProperties::atSw(sw);
    }
    /*!
     * \brief The wetting-phase saturation of a pore body
     *
     * \param pc The capillary pressure \f$\mathrm{p_{c,i}}\f$ at pore \f$i\f$
     * \param params The parameters container
     */
    template<class Scalar,class PcSwDropTables>
    static Scalar sw(Scalar pc, const Params<Scalar>& params)
    {
        if (pc < pcCrit(params))
        {
            return swCrit(params);
        }
        else
        {
            const Scalar sigma = params.surfaceTension();
            const Scalar poreRadius = params.poreInscribedRadius();
            const Scalar poreVolume = (2*poreRadius)*(2*poreRadius)*(2*poreRadius);

            return 1-1/std::cbrt(poreVolume)*4/3*M_PI*(2*sigma/pc)*(2*sigma/pc)*(2*sigma/pc);
        }
        //PcSwTables::TabulatedSw tabulatedSw;
        //return tabulatedSw.TabulatedPcSwProperties::at(pc);
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the wetting phase saturation.
     *
     * \param sw Saturation of the wetting phase \f$\mathrm{S_{w,i}}\f$ at pore \f$i\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar,class PcSwDropTables>
    static Scalar dpc_dsw(Scalar sw, const Params<Scalar>& params)
    {
        PcSwTables::TabulatedPc tabulatedPc;
        Scalar eps = 1e-2;
        if (sw>=0.0+eps)
            return (tabulatedPc.TabulatedPcSwProperties::atSw(sw) - tabulatedPc.TabulatedPcSwProperties::atSw(sw-eps))/eps;
        else
            return (tabulatedPc.TabulatedPcSwProperties::atSw(sw+eps) - tabulatedPc.TabulatedPcSwProperties::atSw(sw))/eps;
    }

    /*!
     * \brief The partial derivative of the wetting phase saturation
     *        w.r.t. the capillary pressure.
     *
     * \param pc The capillary pressure \f$\mathrm{p_{c,i}}\f$ at pore \f$i\f$
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar,class PcSwDropTables>
    static Scalar dsw_dpc(Scalar pc, const Params<Scalar>& params)
    {
        PcSwTables::TabulatedSw tabulatedSw;
        Scalar eps = 1e-2;
        if (pc>=0.0+eps)
            return (tabulatedSw.TabulatedPcSwProperties::at(pc) - tabulatedSw.TabulatedPcSwProperties::at(pc-eps))/eps;
        else
            return (tabulatedSw.TabulatedPcSwProperties::at(pc+eps) - tabulatedSw.TabulatedPcSwProperties::at(pc))/eps;
    }
    /*!
     * \brief The critical saturation of the drop (saturation at max. capillary pressure)
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar swCrit(const Params<Scalar>& params)
    {
        using std::cos;
        const Scalar poreRadius = params.poreInscribedRadius();
        const Scalar throatRadius = params.throatInscribedRadius();
        const Scalar contactAngle = params.criticalContactAngle();
        const Scalar poreVolume = (2*poreRadius)*(2*poreRadius)*(2*poreRadius);
        const Scalar V_cap_crit = M_PI/3*throatRadius*throatRadius*throatRadius*(2+std::cos(M_PI/2-contactAngle))*(1-std::cos(M_PI/2-contactAngle))*(1-std::cos(M_PI/2-contactAngle));// volume of the sphere cap with critical contact angle
        return (poreVolume-V_cap_crit)/poreVolume;
    }
        /*!
     * \brief The critical capillary pressure of the drop (max. capillary pressure)
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     */
    template<class Scalar>
    static Scalar pcCrit(const Params<Scalar>& params)
    {
        using std::cos;
        const Scalar poreRadius = params.poreInscribedRadius();
        const Scalar poreVolume = (2*poreRadius)*(2*poreRadius)*(2*poreRadius);
        const Scalar contactAngle = params.criticalContactAngle();
        const Scalar sigma = params.surfaceTension();
        const Scalar swCrit = swCrit(params);
        const Scalar V_cap = (1-swCrit)*poreVolume;
        const Scalar dropRadius = std::cbrt(3/M_PI *V_cap/(2+std::cos(M_PI/2-contactAngle))*(1-std::cos(M_PI/2-contactAngle))*(1-std::cos(M_PI/2-contactAngle)));

        return 2*sigma/dropRadius;
    }

};

template<class Scalar, class BaseLaw>
class DropLocalRulesRegularization
{
    using ThisType = DropLocalRulesRegularization<Scalar, BaseLaw>;
    using BaseLawParams = typename BaseLaw::template Params<Scalar>;
public:

    template<class MaterialLaw>
    void init(const MaterialLaw* m, const DropParams<Scalar>& p, const std::string& paramGroup = "")
    {
        baseLawParamsPtr_ = &m->basicParams();
    }
    /*!
     * \brief The regularized saturation-capillary pressure curve
     */
    OptionalScalar<Scalar> pc(const Scalar sw) const
    {
        return {};
    }
    /*!
     * \brief The regularized saturation-capillary pressure curve
     */
    OptionalScalar<Scalar> sw(const Scalar pc) const
    {
        return {};
    }
    /*!
     * \brief The regularized partial derivative of the capillary pressure w.r.t. the saturation
     */
    OptionalScalar<Scalar> dpc_dsw(const Scalar sw) const
    {
        return {};
    }
    /*!
     * \brief The regularized partial derivative of the saturation to the capillary pressure
     */
    OptionalScalar<Scalar> dsw_dpc(const Scalar pc) const
    {
        return {};
    }

    template<class SpatialParams, class Element, class SubControlVolume, class ElemSol>
    void updateParams(const SpatialParams& spatialParams,
                      const Element& element,
                      const SubControlVolume& scv,
                      const ElemSol& elemSol)
    {}
private:
    const BaseLawParams* baseLawParamsPtr_;
};
template<Pore::Shape shape, typename Scalar = double>
using DropLocalRulesDefault = SingleShapeTwoPLocalRules<Scalar, DropLocalRules<shape>, DropLocalRulesRegularization<Scalar, DropLocalRules<shape>>>;

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A default configuration without regularization for using the VanGenuchten material law
 */
template<Pore::Shape shape, typename Scalar = double>
using DropLocalRulesNoReg = SingleShapeTwoPLocalRules<Scalar, DropLocalRules<shape>, Dumux::FluidMatrix::NoRegularization>;


} // end namespace Dumux::PoreNetwork::FluidMatrix

#endif
