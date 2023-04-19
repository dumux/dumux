//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <dune/common/indices.hh>
#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>

#include <dune/istl/solvercategory.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/scalarproducts.hh>

#include <dumux/common/initialize.hh>
#include <dumux/linear/scalarproducts.hh>

namespace Dumux::Test {

struct MockCommunication
{
    MockCommunication(Dune::SolverCategory::Category category)
    : category_(category)
    {}

    template<class X, class S>
    void dot(const X& x, const X& y, S& result) const
    {
        Dune::ScalarProduct<X> scalarProduct;
        result = scalarProduct.dot(x, y);
    }

    const Dune::SolverCategory::Category& category() const
    { return category_; }

private:
    Dune::SolverCategory::Category category_;
};

} // end namespace Dumux::Test

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    using Vector1 = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using Vector2 = Dune::BlockVector<Dune::FieldVector<double, 2>>;

    using BlockVector = Dune::MultiTypeBlockVector<Vector1, Vector2>;
    using Comms = std::array<std::shared_ptr<const Test::MockCommunication>, BlockVector::size()>;

    constexpr std::size_t numDofs = 100;
    using namespace Dune::Indices;
    BlockVector b;
    b[_0].resize(numDofs);
    b[_1].resize(2*numDofs);
    b = 1.0;

    using SP = ParallelMultiTypeScalarProduct<BlockVector, Test::MockCommunication>;
    using Cat = Dune::SolverCategory::Category;

    const auto testMixedCategory = [&](Cat cat0, Cat cat1, Cat catCommon)
    {
        Comms comms;
        comms[0] = std::make_shared<Test::MockCommunication>(cat0);
        for (int i = 1; i < BlockVector::size(); ++i)
            comms[i] = std::make_shared<Test::MockCommunication>(cat1);

        SP sp(comms);
        const auto norm2 = sp.dot(b, b);
        const auto norm = sp.norm(b);

        if (sp.category() != catCommon)
            DUNE_THROW(Dune::Exception, "Wrong category");
        if (Dune::FloatCmp::ne<double>(norm2, numDofs*5.0))
            DUNE_THROW(Dune::Exception, "Wrong scalar product " << norm2 << ", ref: " << numDofs*5.0);
        if (Dune::FloatCmp::ne<double>(norm*norm, numDofs*5.0))
            DUNE_THROW(Dune::Exception, "Wrong norm" << norm << ", ref: " << std::sqrt(numDofs*5.0));
    };

    const auto testSameCategory = [&](Cat cat)
    { testMixedCategory(cat, cat, cat); };

    testSameCategory( (Dune::SolverCategory::sequential) );
    testSameCategory( (Dune::SolverCategory::overlapping) );
    testSameCategory( (Dune::SolverCategory::nonoverlapping) );

    testMixedCategory( (Dune::SolverCategory::sequential), (Dune::SolverCategory::overlapping), (Dune::SolverCategory::overlapping) );
    testMixedCategory( (Dune::SolverCategory::sequential), (Dune::SolverCategory::nonoverlapping), (Dune::SolverCategory::overlapping) );

    return 0;
}
