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
 * \ingroup Assembly
 * \brief Flux operator for the TPFA scheme.
 */
#ifndef DUMUX_ASSEMBLY_CC_TPFA_FLUX_OPERATOR_HH
#define DUMUX_ASSEMBLY_CC_TPFA_FLUX_OPERATOR_HH

#include <utility>
#include <concepts>
#include <type_traits>

// forward declarations
namespace Dune {

template<typename K, int rows, int cols>
class FieldMatrix;

} // namespace Dune

namespace Dumux {
namespace Traits {

template<typename T>
struct IsTensor : public std::false_type {};

template<typename K, int rows, int cols>
struct IsTensor<Dune::FieldMatrix<K, rows, cols>>
: public std::true_type
{};

} // namespace Traits


namespace Concepts {

template<typename T>
concept Tensor = Traits::IsTensor<T>::value or std::floating_point<T>;

template<typename T>
concept FluxOperatorTraits = requires(const T& t) {
    { t.tensor() };
    { t.function() } -> std::floating_point;
    Tensor<std::decay_t<decltype(t.tensor())>>;
};

} // namespace Concepts

template<typename T, typename U>
class FluxOperatorTraits
{
public:
    FluxOperatorTraits(T&& t, U&& u)
    : t_(std::move(t))
    , u_(std::move(u))
    {}

    template<typename VolumeVariables> requires(
        std::invocable<const T, const VolumeVariables&>)
    Concepts::Tensor auto tensor(const VolumeVariables& volVars) const
    { return tensorAccess_(volVars); }

    template<typename VolumeVariables> requires(
        std::invocable<const U, const VolumeVariables&>)
    std::floating_point auto function(const VolumeVariables& volVars) const
    { return variableAccess_(volVars); }

private:
    T t_;
    U u_;
};


namespace Concepts {

template<typename T, typename Accessor>
concept VolumeVariablesStorage = requires(const T& t, const Accessor& acc) {
    typename T::VolumeVariables;

    { t.volVars(acc) };
};

template<typename T,
         typename Scvf,
         typename TIJ>
concept TransmissibilityStorage = requires(T& t,
                                           const T& tConst,
                                           const Scvf& scvf,
                                           const TIJ& tij) {
    { t.setTransmissibilities(scvf, tij) };
    { t.getTransmissibilities(scvf) } -> std::convertible_to<TIJ>;
};

} // namespace Concepts











template<typename T,
         typename ElementGeometry,
         typename ElementVariables,
         typename TransmissibilityStorage>
concept CCFluxOperator = requires(const T& t,
                                  const ElementGeometry& eg,
                                  const ElementVariables& ev,
                                  TransmissibilityStorage& ts,
                                  const TransmissibilityStorage& tsConst) {
    typename T::Fluxes;
    typename T::FaceTransmissibilites;
    Concepts::Indexable<T::Fluxes>;
    Concepts::Indexable<T::FaceTransmissibilites>;

    { t.setTransmissibilities(eg, ev, ts) };
    { t(eg, ev, tsConst, std::declval<typename ElementGeometry::SubControlVolumeFace>()) };
};




// flux per unit area of the scvf
template<typename Scalar,
         typename TensorAccess,
         typename VariableAccess>
class CCTpfaFluxOperator
{
public:
    using Flux = std::array<Scalar, 1>;
    using FaceTransmissibilities = std::array<Scalar, 1>;

    CCTpfaFluxOperator(TensorAccess&& tensorAccess,
                       VariableAccess&& variableAccess)
    : tensorAccess_(std::move(tensorAccess))
    , variableAccess_(std::move(variableAccess))
    {}

    template<typename EG, typename EV, typename TS>
    void setTransmissibilities(const EG& eg,
                               const EV& ev,
                               TS& ts) const
    {
        for (const auto& scvf : scvfs(eg))
            ts.setTransmissibilities(
                scvf,
                FaceTransmissibilities{computeTransmissibility_(eg, ev, scvf)}
            );
    }

    template<typename EG,
             typename EV,
             typename TS,
             typename SCVF>
    Scalar operator()(const EG& eg,
                      const EV& ev,
                      const TS& ts,
                      const SCVF& scvf) const
    {
        const FaceTransmissibilities& tij = ts.getTransmissibilities(scvf);

        if (scvf.boundary())
            return {tij[0]*(getVariable_(insideVolVars) - getVariable_(ev.boundaryVolVars(scvf)))};
        else
        {
            const auto& outsideScv = eg.outsideScv(scvf);
            const auto& outsideVolVars = ev.volVars(outsideScv);
            return {tij[0]*(getVariable_(insideVolVars) - getVariable_(outsideVolVars))};
        }
    }

private:
    template<typename EG, typename EV, typename SCVF>
    Scalar computeTransmissibility_(const EG& eg,
                                    const EV& ev,
                                    const SCVF& scvf) const
    {
        const auto& insideScv = eg.insideScv(scvf);
        const auto& insideVolVars = ev.volVars(insideScv);
        const auto ti =  computeTpfaTransmissibility(
            scvf, insideScv, getTensor_(insideVolVars), insideVolVars.extrusionFactor()
        );
        if (scvf.boundary())
            return ti;

        const auto& outsideScv = eg.outsideScv(scvf);
        const auto& outsideVolVars = ev.volVars(outsideScv);
        const auto tj = computeTpfaTransmissibility(
            scvf, outsideScv, getTensor_(outsideVolVars), outsideVolVars.extrusionFactor()
        );
        return ti*tj/(ti - tj);
    }

    template<typename VolVars> requires(
        std::invocable<TensorAccess, VolVars>)
    decltype(auto) getTensor_(const VolVars& volVars) const
    { return tensorAccess_(volVars); }

    template<typename VolVars> requires(
        std::invocable<VariableAccess, VolVars> and
        std::floating_point<std::invoke_result_t<VariableAccess, VolVars>>)
    decltype(auto) getVariable_(const VolVars& volVars) const
    { return variableAccess_(volVars); }

    TensorAccess tensorAccess_;
    VariableAccess variableAccess_;
};

template<typename Scalar,
         typename TensorAccess,
         typename VariableAccess>
auto makeCCTpfaFluxOperator(TensorAccess&& t, VariableAccess&& v)
{
    return CCTpfaFluxOperator<Scalar, TensorAccess, VariableAccess>{
        std::forward<TensorAccess>(t),
        std::forward<VariableAccess>(v)
    };
}


template<typename Scalar>
class CCTpfaDarcysLawFluxOperator
{
    using BaseOperator = decltype(makeBaseOperator_());

public:
    using Flux = std::array<
    using FaceTransmissibilities = typename BaseOperator::FaceTransmissibilities;

    template<typename EG, typename EV, typename TS>
    void setTransmissibilities(const EG& eg,
                               const EV& ev,
                               TS& ts) const
    {
        const auto op = makeBaseOperator_();
        op.setTransmissibilities(eg, ev, ts);
    }

    template<typename EG,
             typename EV,
             typename TS,
             typename SCVF>
    Scalar operator()(const EG& eg,
                      const EV& ev,
                      const TS& ts,
                      const SCVF& scvf) const
    {
        const auto op = makeBaseOperator_();
        return op(eg, ev, ts, scvf);
    }

private:
    static auto makeBaseOperator_() const
    {
        return makeCCTpfaFluxOperator<Scalar>(
            [] (const auto& volVars) { return volVars.permeability(); },
            [] (const auto& volVars) { return volVars.pressure(); }
        );
    }
};
























template<typename UpwindTerm, typename BasicOperator>
class UpwoundFluxOperator
{
public:
    UpwoundFluxOperator(UpwindTerm&& up,
                        BasicOperator&& op)
    : upwindTerm_(std::move(up))
    , operator_(std::move(op))
    {}

    template<typename EG, typename EV, typename TS>
    void setTransmissibilities(const EG& eg,
                               const EV& ev,
                               TS& ts) const
    { operator_.setTransmissibilities(eg, ev, ts); }

    template<typename EG,
             typename EV,
             typename TS,
             typename SCVF>
    Scalar operator()(const EG& eg,
                      const EV& ev,
                      const TS& ts,
                      const SCVF& scvf) const
    {
        const auto flux =
        const FaceTransmissibilities& tij = ts.getTransmissibilities(scvf);

        if (scvf.boundary())
            return {tij[0]*(getVariable_(insideVolVars) - getVariable_(ev.boundaryVolVars(scvf)))};
        else
        {
            const auto& outsideScv = eg.outsideScv(scvf);
            const auto& outsideVolVars = ev.volVars(outsideScv);
            return {tij[0]*(getVariable_(insideVolVars) - getVariable_(outsideVolVars))};
        }
    }

private:
    UpwindTerm upwindTerm_;
    BasicOperator operator_;
};





template<typename Scalar>
class DarcysLawOperator
{
    static constexpr int numOperators = sizeof...(Operators);

    template<int i>
    using OperatorType = std::tuple_element_t<i, std::tuple<Operators...>>;

public:
    using Fluxes = std::array<Scalar, sizeof...(Operators)>;
    using FaceTransmissibilities = typename OperatorType<0>::FaceTransmissibilities;

    SharedTransmissibilitiesOperator(Operators&&... ops)
    : operators_(std::forward<Operators>(ops)...)
    {}

    template<typename EG, typename EV, typename TS>
    void setTransmissibilities(const EG& eg,
                               const EV& ev,
                               TS& ts) const
    {
        std::get<0>(operators_).setTransmissibilities(eg, ev, ts);
    }

    template<typename EG,
             typename EV,
             typename TS,
             typename SCVF>
    Fluxes operator()(const EG& eg,
                      const EV& ev,
                      const TS& ts,
                      const SCVF& scvf) const
    {
        Fluxes result;
        addFluxes_<0>(result, eg, ev, ts, scvf);
        return result;
    }

private:
    template<int curOperatorIdx,
             typename EG,
             typename EV,
             typename TS,
             typename SCVF>
    void addFluxes_(Fluxes& fluxes,
                    const EG& eg,
                    const EV& ev,
                    const TS& ts,
                    const SCVF& scvf) const
    {
        fluxes[curOperatorIdx] = std::get<curOperatorIdx>(operators_)(eg, ev, ts, scvf);
        if constexpr (curOperatorIdx < numOperators - 1)
            addFluxes_<curOperatorIdx+1>(fluxes, eg, ev, ts, scvf);
    }

private:
    std::tuple<Operators...> operators_;
};





template<typename Scalar>
class CCTpfaDarcysLawOperator
{
public:
    template<typename EG, typename EV, typename SCVF>
    std::floating_point auto operator()(const EG& eg,
                                        const EV& ev,
                                        const SCVF& scvf) const
    {
        CCTpfaFluxOperator{
            [] (auto& volVars) { return volVars.permeability(); },
            [] (auto& volVars) { return volVars.pressure(0); }
        }(eg, ev, scvf);
    }
};

template<int fluxDim, typename Operator>
class OperatorMap
{
public:
    OperatorMap(std::array<int, fluxDim>&& indices,
                Operator&& op)
    : indices_(indices)
    , op_(op)
    {}

private:

};

struct StorageTermOneP
{
    template<typename VolumeVariables>
    auto operator()(const VolumeVariables& volVars) const
    { return volVars.porosity()*volVars.density(); }
};

struct FluxTermOneP
{
    template<typename VolumeVariables>
    auto operator()(const VolumeVariables& volVars) const
    { return volVars.porosity()*volVars.density(); }
};

template<typename Scalar,
         typename FluidSystem>
class OnePModel
{
    using VolumeVariables = void;

public:
    struct Indices
    {
        static constexpr int pressureIdx = 0;
        static constexpr int conti0EqIdx = 0;
    };

    using StorageOperators = void;
    using FluxOperators = std::tuple<
        CCTpfaDarcysLawOperator<Scalar>
    >;


};

} // namespace Dumux

#endif
