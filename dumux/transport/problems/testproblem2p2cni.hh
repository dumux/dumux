//$Id$
#ifndef TESTPROBLEM2P2CNI_HH
#define TESTPROBLEM2P2CNI_HH

#include <iostream>
#include <dumux/transport/transportproblem2p2cni.hh>

namespace Dune
{
//! Base class for the definition of nonisothermal 2p2c problems
/** This base class defines all boundary and initial functions which are needed
 * for a decoupled nonisothermal 2p2c computation.
 */
template<class G, class RT>
class TestProblem2p2cni: public TransportProblem2p2cni<G,RT>
{

    typedef typename G::ctype DT;
    enum {n=G::dimension, m=1, blocksize=2*G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;

public:

    virtual BoundaryConditions2p2c::Flags cbctype (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return BoundaryConditions2p2c::concentration;
    }

    virtual BoundaryConditions2p2c::Flags ictype (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        if (x[0] > 60 && x[0] < 120 && x[1] > 120 && x[1] < 180) return BoundaryConditions2p2c::saturation;
        return BoundaryConditions2p2c::concentration;
    }

    virtual BoundaryConditions::Flags pbctype (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        if (x[0] < 1e-6 || x[0] > 300 - 1e-6)
            return BoundaryConditions::dirichlet;
        return BoundaryConditions::neumann;
    }

    virtual RT gZ (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return 1.;
    }

    virtual RT gS (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return 0.1;
    }

    virtual RT gPress (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        if (x[0] < 1e-6) return 2e5;
        if (x[0] > 300 - 1e-6) return 1e5;
    }

    virtual RT gT (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return 290.15;
    }

    virtual FieldVector<RT,3> J (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,3> J_(0);
        if (x[0]<1e-6)
        {
            J_[2] = 10;
            J_[0] = 1;
        }
        return J_;
    }

    virtual FieldVector<RT,2> q (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,2> q_(0);
        return q_;
    }

    virtual RT qh (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        if (x[0] > 60 && x[0] < 120 && x[1] > 120 && x[1] < 180) return 100.;
        return 0.;
    }

    virtual RT S0 (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return 0.1;
    }

    virtual RT Z1_0 (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return 1.;
    }

    virtual RT T_0 (const FieldVector<DT,n>& x, const Entity& e, const FieldVector<DT,n>& xi) const
    {
        return 283.15;
    }

    const FieldVector<DT,n>& gravity()
    {
        return this->gravity_;
    }


    //! Constructor
    /**
     *
     */
    TestProblem2p2cni(Dune::VariableClass2p2cni<G, RT>& var, Liquid_GL& liq, Gas_GL& gas, Matrix2p<G, RT>& s, TwoPhaseRelations<G, RT>& law, const bool cap = false)
        :TransportProblem2p2cni<G,RT>(var, liq, gas, s, law, cap)
    {
        this->gravity_ = 0.;
    }

};

}
#endif
