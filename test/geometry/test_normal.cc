// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for normal computation
 */
#include <config.h>

#include <random>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/geometry/normal.hh>

template<class P>
void testOrthogonality(const P& normal, const P& vector, double scale, const std::string& caseName = "")
{
    const auto eps = Dune::FloatCmp::DefaultEpsilon<double, Dune::FloatCmp::absolute>::value()*scale*100;
    if (Dune::FloatCmp::ne<double, Dune::FloatCmp::absolute>(normal*vector, 0.0, eps))
        DUNE_THROW(Dune::Exception, "Normal not orthogonal, n: " << normal << ", p: " << vector
                    << ", dim: " << P::size() << ", sp: " << normal*vector
                    << ", scale: " << scale << ", epsilon: " << eps << ", case: " << caseName);
}

template<int dim, std::size_t repetitions = 100000>
void testNormal()
{
    std::default_random_engine gen(std::random_device{}());
    std::uniform_real_distribution<double> rnd(-1.0, 1.0), rndIdx(0.0, dim-1);
    Dune::FieldVector<double, dim> p(0.0);

    for (double scale : {1.0, 1e3, 1e-12, 1e12})
    {
        for (std::size_t i = 0; i < repetitions; ++i)
        {
            std::generate(p.begin(), p.end(), [&]{ return rnd(gen)*scale; });
            const auto n = Dumux::normal(p);
            testOrthogonality(n, p, scale, "rnd");
        }

        // one coordinate != 0
        for (std::size_t i = 0; i < repetitions; ++i)
        {
            p = 0.0;
            const auto randomIndex = (std::size_t) std::round(rndIdx(gen));
            p[randomIndex] = scale;
            testOrthogonality(Dumux::normal(p), p, scale, "rnd-one-nonzero");
            p[randomIndex] = scale*1e-20;
            testOrthogonality(Dumux::normal(p), p, scale, "rnd-one-small");
        }

        // two coordinates != 0
        for (std::size_t i = 0; i < repetitions; ++i)
        {
            p = 0.0;
            const auto randomIndex0 = (std::size_t) std::round(rndIdx(gen));
            const auto randomIndex1 = (std::size_t) std::round(rndIdx(gen));
            p[randomIndex0] = scale;
            p[randomIndex1] = scale;
            testOrthogonality(Dumux::normal(p), p, scale, "rnd-two-nonzero");
            p[randomIndex0] = scale*1e-20;
            p[randomIndex1] = scale*1e-20;
            testOrthogonality(Dumux::normal(p), p, scale, "rnd-two-small");
        }
    }
}

int main(int argc, char** argv)
{
    testNormal<2>();
    testNormal<3>();
    testNormal<4>();
    testNormal<10>();

    return 0;
}
