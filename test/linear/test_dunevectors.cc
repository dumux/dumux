//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>

#include <dumux/linear/dunevectors.hh>

int main(int argc, char* argv[])
{
    // scalar type ///////
    static_assert(std::is_same_v<double, typename Dumux::Detail::NativeDuneVectorType<double>::type>);

    // vector types ///////
    using Block1 = Dune::FieldVector<double, 1>;
    using Block2 = Dune::FieldVector<double, 2>;
    using Block3 = double;

    using SwitchableBlock1 = Dumux::SwitchablePrimaryVariables<Block1, int>;
    using SwitchableBlock2 = Dumux::SwitchablePrimaryVariables<Block2, int>;

    using Vector1 = Dune::BlockVector<Block1>;
    using Vector2 = Dune::BlockVector<Block2>;
    using Vector3 = Dune::BlockVector<Block3>;

    using Vector1S = Dune::BlockVector<SwitchableBlock1>;
    using Vector2S = Dune::BlockVector<SwitchableBlock2>;

    // strip switchable
    static_assert(std::is_same_v<Vector1, typename Dumux::Detail::NativeDuneVectorType<Vector1S>::type>);
    static_assert(std::is_same_v<Vector2, typename Dumux::Detail::NativeDuneVectorType<Vector2S>::type>);

    // identity maps
    static_assert(std::is_same_v<Vector1, typename Dumux::Detail::NativeDuneVectorType<Vector1>::type>);
    static_assert(std::is_same_v<Vector2, typename Dumux::Detail::NativeDuneVectorType<Vector2>::type>);
    static_assert(std::is_same_v<Vector3, typename Dumux::Detail::NativeDuneVectorType<Vector3>::type>);

    // multi-type vector types ///////
    using BlockVector11 = Dune::MultiTypeBlockVector<Vector1, Vector1>;
    using BlockVector12 = Dune::MultiTypeBlockVector<Vector1, Vector2>;
    using BlockVector13 = Dune::MultiTypeBlockVector<Vector1, Vector3>;
    using BlockVector23 = Dune::MultiTypeBlockVector<Vector2, Vector3>;
    using BlockVector1S1S = Dune::MultiTypeBlockVector<Vector1S, Vector1S>;
    using BlockVector1S2 = Dune::MultiTypeBlockVector<Vector1S, Vector2>;
    using BlockVector1S2S = Dune::MultiTypeBlockVector<Vector1S, Vector2S>;
    using BlockVector2S3 = Dune::MultiTypeBlockVector<Vector2S, Vector3>;

    // strip switchable
    static_assert(std::is_same_v<BlockVector11, typename Dumux::Detail::NativeDuneVectorType<BlockVector1S1S>::type>);
    static_assert(std::is_same_v<BlockVector12, typename Dumux::Detail::NativeDuneVectorType<BlockVector1S2>::type>);
    static_assert(std::is_same_v<BlockVector12, typename Dumux::Detail::NativeDuneVectorType<BlockVector1S2S>::type>);
    static_assert(std::is_same_v<BlockVector23, typename Dumux::Detail::NativeDuneVectorType<BlockVector2S3>::type>);

    // identity maps
    static_assert(std::is_same_v<BlockVector11, typename Dumux::Detail::NativeDuneVectorType<BlockVector11>::type>);
    static_assert(std::is_same_v<BlockVector12, typename Dumux::Detail::NativeDuneVectorType<BlockVector12>::type>);
    static_assert(std::is_same_v<BlockVector13, typename Dumux::Detail::NativeDuneVectorType<BlockVector13>::type>);

    return 0;
}
