// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Coupling reconstruction for 1d-3d mixed-dimension method with distributed sources
 */
#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLING_RECONSTRUCTION_1D3D_KERNEL_ROOTSOIL_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLING_RECONSTRUCTION_1D3D_KERNEL_ROOTSOIL_HH

#include <cmath>
#include <algorithm>
#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/integrate.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>

#include <dumux/nonlinear/findscalarroot.hh>

namespace Dumux::Detail::RootSoil {

/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Kirchhoff transformation for the Richards equation
 */
class KirchhoffTransformation
{
    double alpha, n, m, l;
public:
    KirchhoffTransformation()
    {
        alpha = getParam<double>("Soil.SpatialParams.VanGenuchtenAlpha");
        n = getParam<double>("Soil.SpatialParams.VanGenuchtenN");
        l = getParam<double>("Soil.SpatialParams.VanGenuchtenL", 0.5);
        m = 1.0 - 1.0/n;
    }

    inline double kr(const double pw) const
    {
        using std::pow;
        const auto pc = 1.0e5-pw;
        const auto swe = pow(pow(alpha*pc, n) + 1.0, -m);
        const auto r = 1.0 - pow(1.0 - pow(swe, 1.0/m), m);
        return pow(swe, l)*r*r;
    }

    //! Kirchhoff transform
    inline double computeTransformed(const double p) const
    { return integrateScalarFunction([this](const auto& p){ return kr(p); }, 0.0, p, 1e-15); }

    //! Inverse Kirchhoff transform
    inline double computeInverseTransformed(const double psi) const
    {
        const auto residual = [this, psi](const auto& p){ return psi - computeTransformed(p); };
        return findScalarRootBrent(-2.6e7, 1.0e5, residual, 1e-15, 200000);
    }
};

/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Kirchhoff transformation for the Richards equation with values cached in a lookup table
 */
class KirchhoffTransformationCached
{
    double pMin, pMax;
    std::size_t samples;
    KirchhoffTransformation transformation_;
    std::vector<double> p_, psi_;
public:
    KirchhoffTransformationCached()
    {
        pMin = getParam<double>("MixedDimension.Reconstruction.PressureMin", -2.6e7);
        pMax = getParam<double>("MixedDimension.Reconstruction.PressureMax", 1.0e5);
        samples = getParam<std::size_t>("MixedDimension.Reconstruction.Samples", 10000);

        p_.resize(samples);
        psi_.resize(samples);
        const double pStep = (pMax-pMin)/samples;
        for (std::size_t i = 0; i <= samples; ++i)
        {
            const auto pValue = pMin + i*pStep;
            p_[i] = pValue;
            psi_[i] = transformation_.computeTransformed(pValue);
        }

        std::cout << "[coupling] Presampled Kirchhoff transformation with "
                  << samples << " samples (step size: " << pStep << " Pa)."
                  << std::endl;
    }

    //! Kirchhoff transform
    inline double computeTransformed(const double p) const
    { return interpolate<InterpolationPolicy::LinearTable>(p, p_, psi_); }

    //! Inverse Kirchhoff transform
    inline double computeInverseTransformed(const double psi) const
    { return interpolate<InterpolationPolicy::LinearTable>(psi, psi_, p_); }
};

/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Caching the source computation to speed up Jacobian approximation during assembly
 */
struct SourceCache
{
    SourceCache()
    : q_{}
    {
        p0_.fill(1e20);
        pRoot_.fill(1e20);
    }

    template<class F>
    double operator() (const double p0, const double pRoot, const F& f)
    {
        for (std::size_t i = 0; i < size; ++i)
            if (std::abs(p0-p0_[i]) < 1e-14 && std::abs(pRoot-pRoot_[i]) < 1e-14)
                return q_[i];

        const auto q = f(p0, pRoot);
        p0_[writePos] = p0;
        pRoot_[writePos] = pRoot;
        q_[writePos] = q;
        writePos = (writePos >= size-1) ? 0 : writePos+1;
        return q;
    }

private:
    static constexpr std::size_t size = 3; // original and each variable deflected
    std::array<double, size> q_, p0_, pRoot_;
    std::size_t writePos = 0;
};

} // end namespace Dumux::Detail::RootSoil

namespace Dumux::RootSoil {

/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Coupling reconstruction as derived in Koch et al (2022), https://doi.org/10.1016/j.jcp.2021.110823
 */
template<class PcKrSwCurve>
class CouplingReconstruction
{
    using CacheStorage = std::vector<std::vector<std::pair<std::size_t, Detail::RootSoil::SourceCache>>>;
public:
    CouplingReconstruction()
    {
        K_ = getParam<double>("Soil.SpatialParams.Permeability");
        porosity_ = getParam<double>("Soil.SpatialParams.Porosity");
        // relative tolerance in term of the unknown approximation (psBar)
        nonlinearSolverTolerance_ = getParam<double>("MixedDimension.Reconstruction.RelativeTolerance", 1e-4);

        std::cout << "\n" << "[coupling] Test Kirchhoff transform:\n";
        const auto& tf = transformation_;
        for (const auto p : {0.9e5, 0.5e5, 0.1e5, 0.0, -1.0e5, -5.0e5, -1.5e6})
        {
            const auto psi = tf.computeTransformed(p);
            const auto pNum = tf.computeInverseTransformed(psi);
            std::cout << "pw: " << p << ", psi: " << psi << ", pw (inverse): " << pNum << '\n';
            if (std::abs(p-pNum) > 100)
                std::cerr << "Transformation is not precise enough. Maybe decrease the tolerances!" << std::endl;
        }
        std::cout << std::endl;
    }

    double computePsBar(const double p0, const double pRoot,
                        const double kernelRadius, const double rootRadius,
                        const double Kr, const double density, const double viscosity,
                        const std::size_t sourceId,
                        const double delta) const
    {
        const auto& tf = transformation_;
        const auto psi0 = tf.computeTransformed(p0);
        const auto sourceFactor = (std::log(kernelRadius/rootRadius) - 0.5 + 0.5*delta*delta/kernelRadius/kernelRadius)*viscosity/(2.0*M_PI*K_);

        const auto residual = [&](const double psBar){
            const auto sourceTerm = sourceFactor*2.0*M_PI*rootRadius*Kr*(psBar - pRoot);
            return psi0 - tf.computeTransformed(psBar) - sourceTerm;
        };

        return findScalarRootBrent(-2.5e7, 1.0e5, residual, nonlinearSolverTolerance_, 20000);
    }

    double computeSource(const double p0, const double pRoot,
                         const double kernelRadius, const double rootRadius,
                         const double Kr, const double density, const double viscosity,
                         const std::size_t sourceId,
                         const double delta) const
    {
        const auto psBar = computePsBar(p0, pRoot, kernelRadius, rootRadius, Kr, density, viscosity, sourceId, delta);
        return -2.0*M_PI*rootRadius*Kr*(psBar - pRoot)*density;
    }

    // try to be smart and use a cache to reduce nonlinear solves
    double computeSource(std::size_t bEIdx, std::size_t lEIdx,
                         const double p0, const double pRoot,
                         const double kernelRadius, const double rootRadius,
                         const double Kr, const double density, const double viscosity,
                         const std::size_t sourceId,
                         const double delta) const
    {
        if (!sourceCache_)
            return computeSource(p0, pRoot, kernelRadius, rootRadius, Kr, density, viscosity, sourceId, delta);

        auto& cacheLine = (*sourceCache_)[lEIdx];
        auto bIt = std::find_if(cacheLine.begin(), cacheLine.end(), [&](const auto& p){ return p.first == bEIdx; });
        if (bIt == cacheLine.end())
            DUNE_THROW(Dune::InvalidStateException, "Source cache does not exist for index pair (l,b): " << lEIdx << "," << bEIdx);

        return bIt->second(p0, pRoot, [&](const auto p0_, const auto pRoot_)
        {
            const auto psBar = computePsBar(p0_, pRoot_, kernelRadius, rootRadius, Kr, density, viscosity, sourceId, delta);
            return -2.0*M_PI*rootRadius*Kr*(psBar - pRoot_)*density;
        });
    }

    const  Dumux::Detail::RootSoil::KirchhoffTransformationCached& transformation() const
    { return transformation_; }

    template<class CouplingManager>
    void enableCache(const CouplingManager& couplingManager)
    {
        // create a cache vector for each low dim element
        // the vector shall contain one element for each coupled bulk element (where p0 is measured)
        const auto& lowDimGG = couplingManager.problem(CouplingManager::lowDimIdx).gridGeometry();
        sourceCache_ = std::make_unique<CacheStorage>(lowDimGG.numDofs());
        // lets use the stencil to find these coupled elements (probably only works for tpfa)
        for (std::size_t lowDimElementIdx = 0; lowDimElementIdx < sourceCache_->size(); ++lowDimElementIdx)
        {
            const auto& couplingStencil = couplingManager.couplingStencils(CouplingManager::lowDimIdx).at(lowDimElementIdx);
            (*sourceCache_)[lowDimElementIdx].reserve(couplingStencil.size());
            for (std::size_t bulkElementIdx : couplingStencil)
                (*sourceCache_)[lowDimElementIdx].emplace_back(bulkElementIdx, Detail::RootSoil::SourceCache{});
        }
    }
private:
    std::unique_ptr<CacheStorage> sourceCache_;

    Detail::RootSoil::KirchhoffTransformationCached transformation_;
    double K_, porosity_;
    double nonlinearSolverTolerance_;
};

} // end namespace Dumux::RootSoil

#endif
