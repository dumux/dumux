#ifndef DUNE_SEQCOUPBASEPROBLEM_HH
#define DUNE_SEQCOUPBASEPROBLEM_HH

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/auxiliary/basicdomain.hh>

#include "seqcoup_2pni_problem.hh"
#include "seqcoup_2p2cni_problem.hh"
#include "seqcoup_base_problem.hh"
/**
 * @file
 * @brief  Definition of a problem, where a 2p2cni model is sequentially coupled to a 2pni model
 * @author Bernd Flemisch, Melanie Darcis
 */

namespace Dune
{
//! class that defines the parameters of an air waterair under a low permeable layer
/*! Problem definition of an air injection under a low permeable layer. Air enters the domain
 * at the right boundary and migrates upwards.
 * Problem was set up using the rect2d.dgf grid.
 *
 *    Template parameters are:
 *
 *    - ScalarT  Floating point type used for scalars
 */
template<class ScalarT>
class SeqCoupBaseProblem : public BasicDomain<Dune::YaspGrid<2>,
                                              ScalarT>
{
    typedef YaspGrid<2>                   Grid;
    typedef BasicDomain<Grid, ScalarT>          ParentType;
    typedef SeqCoup2PNIProblem<ScalarT>         Problem2PNI;
    typedef SeqCoup2P2CNIProblem<ScalarT>       Problem2P2CNI;
    typedef GridPtr<Grid>                 GridPointer;

public:
    // the domain traits of the domain
    typedef typename ParentType::DomainTraits   DomainTraits;
    // copy some types from the traits for convenience
    typedef typename DomainTraits::Scalar                     Scalar;
    typedef typename Dune::BlockVector<Dune::FieldVector<double,3> > TwoPNISolution;
    typedef typename Problem2PNI::Model TwoPNIModel;

private:
    // some constants from the traits for convenience
    enum {
        // Grid and world dimension
        dim  = DomainTraits::dim,
        dimWorld = DomainTraits::dimWorld,
    };
public:
    SeqCoupBaseProblem(Grid *grid, Scalar tEnd2pni, Scalar dt2pni, Scalar tEnd2p2cni, Scalar dt2p2cni)
        : ParentType(grid),
          twoPNISolution_(this->numVertices()),
          problem2pni_(grid, tEnd2pni, dt2pni),
          problem2p2cni_(grid, tEnd2p2cni, dt2p2cni)
    {

    }


    bool simulate()
    {

        problem2pni_.simulate();
        //Lösungsvector raussschreiben
        // twoPNISolution_ =
        twoPNISolution_ = (*problem2pni_.model().currentSolution());
        problem2p2cni_.sequentialCoupling(twoPNISolution_);
        problem2p2cni_.simulate();
        return true;
    };

    void deserialize(double t)
    {
        problem2pni_.deserialize(t);
    };

private:
    TwoPNISolution twoPNISolution_;

    Problem2PNI problem2pni_;
    Problem2P2CNI problem2p2cni_;
    bool wasRestarted_;
    bool sequentialCoupling_;
};
} //end namespace

#endif
