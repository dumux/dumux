// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Some tools for random number generation
 */
#ifndef DUMUX_COMMON_RANDOM_HH
#define DUMUX_COMMON_RANDOM_HH

#include <random>
#include <type_traits>
#include <cstdint>

namespace Dumux {

/*!
 * \brief A simple uniform distribution
 * based on a biased uniform number generator
 * \note Use this if you need a fast library implementation independent generator
 *       without strict requirements about the bias
 * \note We try to stay close to https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
 */
template<class Scalar = double>
class SimpleUniformDistribution
{
    struct Parameters
    {
        Parameters(Scalar a, Scalar b)
        : a_(a), b_(b) {}

        Scalar a() const { return a_; }
        Scalar b() const { return b_; }
    private:
        Scalar a_, b_;
    };
public:
    using param_type = Parameters;
    using result_type = Scalar;

    explicit SimpleUniformDistribution(Scalar min, Scalar max = 1.0)
    : offset_(min)
    , size_(max-min)
    {}

    explicit SimpleUniformDistribution(const Parameters& p)
    : SimpleUniformDistribution(p.a(), p.b())
    {}

    SimpleUniformDistribution()
    : SimpleUniformDistribution(0.0)
    {}

    void param(const Parameters& p)
    {
        offset_ = p.a();
        size_ = p.b()-p.a();
    }

    Parameters param() const
    { return { offset_, offset_+size_ }; }

    Scalar a() const { return offset_; }
    Scalar b() const { return offset_ + size_; }

    template<class Generator,
             typename std::enable_if_t<std::is_same_v<typename Generator::result_type, std::uint_fast32_t>, int> = 0>
    Scalar operator()(Generator& gen)
    { return offset_ + size_*(0x1.0p-32 * gen()); }

private:
    Scalar offset_;
    Scalar size_;
};

/*!
 * \brief A simple normal distribution
 * based on a biased uniform number generator and the Box-Mueller transform
 * \note Use this if you need a fast library implementation independent generator
 *       without strict requirements about the bias
 * \note We try to stay close to https://en.cppreference.com/w/cpp/numeric/random/normal_distribution
 */
template<class Scalar = double>
class SimpleNormalDistribution
{
    struct Parameters
    {
        Parameters(Scalar m, Scalar s)
        : m_(m), s_(s) {}

        Scalar m() const { return m_; }
        Scalar s() const { return s_; }
    private:
        Scalar m_, s_;
    };
public:
    using param_type = Parameters;
    using result_type = Scalar;

    explicit SimpleNormalDistribution(Scalar mean, Scalar stddev = 1.0)
    : mean_(mean)
    , stddev_(stddev)
    {}

    explicit SimpleNormalDistribution(const Parameters& p)
    : SimpleNormalDistribution(p.m(), p.s())
    {}

    SimpleNormalDistribution()
    : SimpleNormalDistribution(0.0)
    {}

    void param(const Parameters& p)
    {
        mean_ = p.m();
        stddev_ = p.s();
    }

    Parameters param() const
    { return { mean_, stddev_ }; }

    Scalar m() const { return mean_; }
    Scalar s() const { return stddev_; }

    template<class Generator>
    Scalar operator()(Generator& gen)
    {
        if (isCached_)
        {
            isCached_ = false;
            return cachedValue_;
        }

        // Box-Mueller transform (https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform)
        const auto [u1, u2] = generateUniformPair_(gen);

        using std::sqrt; using std::log;
        using std::cos; using std::sin;
        constexpr Scalar twoPi = 2.0 * M_PI;

        const Scalar magnitude = stddev_ * sqrt(-2.0 * log(u1));
        const Scalar z0  = magnitude * cos(twoPi * u2) + mean_;
        const Scalar z1  = magnitude * sin(twoPi * u2) + mean_;
        cachedValue_ = z0;
        isCached_ = true;
        return z1;
    }

private:
    template<class Generator,
             typename std::enable_if_t<std::is_same_v<typename Generator::result_type, std::uint_fast32_t>, int> = 0>
    auto generateUniformPair_(Generator& gen)
    {
        // biased uniform number generator (0,1)
        // https://www.pcg-random.org/posts/bounded-rands.html
        constexpr Scalar eps = std::numeric_limits<Scalar>::epsilon();
        Scalar u1 = 0.0, u2 = 0.0;
        do {
            u1 = 0x1.0p-32 * gen();
            u2 = 0x1.0p-32 * gen();
        } while (u1 <= eps);
        return std::make_pair(u1, u2);
    }

    Scalar mean_;
    Scalar stddev_;
    bool isCached_ = false;
    Scalar cachedValue_ = {};
};

/*!
 * \brief A simple log-normal distribution
 * \note Use this if you need a fast library implementation independent generator
 *       without strict requirements about the bias
 * \note We try to stay close to https://en.cppreference.com/w/cpp/numeric/random/lognormal_distribution
 */
template<class Scalar = double>
class SimpleLogNormalDistribution
{
    using Parameters = typename SimpleNormalDistribution<Scalar>::param_type;
public:
    using param_type = Parameters;
    using result_type = Scalar;

    explicit SimpleLogNormalDistribution(Scalar mean, Scalar stddev = 1.0)
    : normal_(mean, stddev)
    {}

    explicit SimpleLogNormalDistribution(const Parameters& p)
    : SimpleLogNormalDistribution(p.mean, p.stddev)
    {}

    SimpleLogNormalDistribution()
    : SimpleLogNormalDistribution(0.0)
    {}

    void param(const Parameters& p)
    { normal_.param(p); }

    Parameters param() const
    { return normal_.param(); }

    Scalar m() const { return normal_.m(); }
    Scalar s() const { return normal_.s(); }

    template<class Generator>
    Scalar operator()(Generator& gen)
    {
        using std::exp;
        return exp(normal_(gen));
    }

private:
    SimpleNormalDistribution<Scalar> normal_;
};

} // end namespace Dumux

#endif
