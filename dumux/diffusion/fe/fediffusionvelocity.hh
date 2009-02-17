#ifndef DUNE_FEDIFFUSION_HH
#define DUNE_FEDIFFUSION_HH

namespace Dune
{

template<class G, class RT, class VC, class Problem = FractionalFlowProblem<G, RT, VC>, class LocalStiffnessType = GroundwaterEquationLocalStiffness<G,RT,Problem> >
void FEDiffusion<G, RT, VC, Problem, LocalStiffnessType>::totalVelocity(VelType& velocity, const RT t=0) const
{
	return;
}

}
#endif
