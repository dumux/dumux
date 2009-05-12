// $Id$

#ifndef DUMUX_NEW_1P2C_BOX_MODEL_HH
#define DUMUX_NEW_1P2C_BOX_MODEL_HH

#include <dumux/new_models/1p2c/1p2cboxjacobian.hh>

namespace Dune
{
/**
 * \brief Isothermal one phase two component model with pressure and
 *        molefraction as primary unknowns.
 *
 * This implements an isothermal one phase two component model
 * with pressure and molefraction as primary unknowns.
 */
template<class ProblemT, class OnePTwoCTraitsT = OnePTwoCTraits<typename ProblemT::DomainTraits::Scalar> >
class OnePTwoCBoxModel
    : public BoxScheme<OnePTwoCBoxModel<ProblemT, OnePTwoCTraitsT>, // Implementation of the box scheme

                       // The Traits for the BOX method
                       P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                   typename ProblemT::DomainTraits::Grid,
                                   OnePTwoCTraitsT::numEq>,

                       // The actual problem we would like to solve
                       ProblemT,
                       // The local jacobian operator
                       OnePTwoCBoxJacobian<ProblemT,
                                           P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                       typename ProblemT::DomainTraits::Grid,
                                                       OnePTwoCTraitsT::numEq>,
                                           OnePTwoCTraitsT,
                                           OnePTwoCElementData<OnePTwoCTraitsT,ProblemT>,
                                           OnePTwoCVertexData<OnePTwoCTraitsT,ProblemT> > >
{
    typedef typename ProblemT::DomainTraits::Grid       Grid;
    typedef typename ProblemT::DomainTraits::Scalar     Scalar;
    typedef OnePTwoCBoxModel<ProblemT,OnePTwoCTraitsT>  ThisType;

public:
    typedef OnePTwoCTraitsT                                  OnePTwoCTraits;
    typedef P1BoxTraits<Scalar, Grid, OnePTwoCTraits::numEq> BoxTraits;
    typedef OnePTwoCElementData<OnePTwoCTraitsT,ProblemT> 	 ElementData;
    typedef OnePTwoCVertexData<OnePTwoCTraitsT,ProblemT>	 VertexData;


private:
    typedef OnePTwoCBoxJacobian<ProblemT, BoxTraits, OnePTwoCTraits, ElementData, VertexData>  OnePTwoCLocalJacobian;
    typedef BoxScheme<ThisType,
                      BoxTraits,
                      ProblemT,
                      OnePTwoCLocalJacobian>        ParentType;

    typedef typename ProblemT::DomainTraits              DomTraits;
    typedef typename DomTraits::Element                  Element;
    typedef typename DomTraits::Vertex                   Vertex;
    typedef typename DomTraits::ElementIterator          ElementIterator;
    typedef typename DomTraits::LocalPosition            LocalPosition;
    typedef typename DomTraits::GlobalPosition           GlobalPosition;

    enum {
        dim              = DomTraits::dim,
        dimWorld         = DomTraits::dimWorld
    };

public:
    typedef NewNewtonMethod<ThisType> NewtonMethod;

    OnePTwoCBoxModel(ProblemT &prob)
        : ParentType(prob, onePTwoCLocalJacobian_),
          onePTwoCLocalJacobian_(prob)
    {
        Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
    }


    /*!
     * \brief Called by the update() method if applying the newton
     *         method was unsuccessful.
     */
    void updateFailedTry()
    {
        ParentType::updateFailedTry();

    };

    /*!
     * \brief Called by the BoxScheme's update method.
     */
    void updateSuccessful()
    {
        ParentType::updateSuccessful();

    }


    /*!
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer)
    {
        onePTwoCLocalJacobian_.addVtkFields(writer, this->currentSolution());
    }




private:
    // calculates the jacobian matrix at a given position
    OnePTwoCLocalJacobian  onePTwoCLocalJacobian_;
};
}

#endif

