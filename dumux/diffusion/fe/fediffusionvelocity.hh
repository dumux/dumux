#ifndef DUNE_FEDIFFUSION_HH
#define DUNE_FEDIFFUSION_HH

namespace Dune
{

template<class GridView, class Scalar, class VC, class Problem, class LocalStiffnessType>
void FEPressure2P<GridView, Scalar, VC, Problem, LocalStiffnessType>::calculateVelocity(const Scalar t=0) const
{
    DUNE_THROW(Dune::NotImplemented, "velocities only implemented for finite volume and mimetic finite differences discretisations");
    return;
}

}
#endif
