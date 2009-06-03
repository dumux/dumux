// $Id$

#ifndef DUNE_CAPILLARYDIFFUSION_HH
#define DUNE_CAPILLARYDIFFUSION_HH

#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/transport/transportproblem.hh"

//! \ingroup transport
//! \defgroup diffPart Diffusive transport
/**
 * @file
 * @brief  Base class for defining the diffusive part of an advection-diffusion equation
 * @author Bernd Flemisch, Markus Wolff
 */
namespace Dune
{
/*!\ingroup diffPart
 * @brief  Base class for defining the diffusive part of an advection-diffusion equation
 */
template<class GridView, class Scalar, class VC,
        class Problem = TransportProblem<GridView, Scalar, VC> >
class CapillaryDiffusion: public DiffusivePart<GridView, Scalar>
{
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
typedef    typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    virtual FieldVector operator() (const Element& element, const int indexInInside, Scalar satI, Scalar satJ, const FieldVector& pcGradient) const
    {
        // cell geometry type
        GeometryType gt = element.geometry().type();

        // cell center in reference element
        const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0,0);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = element.geometry().global(localPos);

        // get absolute permeability of cell
        FieldMatrix permeability(soil_.K(globalPos,element,localPos));

        IntersectionIterator isItEnd = element.ilevelend();
        IntersectionIterator isIt = element.ilevelbegin();
        for (; isIt != isItEnd; ++isIt)
        {
            if(isIt->indexInInside() == indexInInside)
            break;
        }
        int globalIdxI = problem_.variables().indexTransport(element);

        // get geometry type of face
        GeometryType faceGT = isIt->geometryInInside().type();

        //get lambda_bar = lambda_n*f_w
        Scalar mobBarI, mobBarJ;

	if (preComput_)
            mobBarI=problem_.variables().mobilityWetting()[globalIdxI]*problem_.variables().fracFlowFuncNonWetting()[globalIdxI];
	else
	    mobBarI=problem_.materialLaw().mobN(1-satI,globalPos,element,localPos)*problem_.materialLaw().fractionalW(satI,globalPos,element,localPos);

        if (isIt->neighbor())
        {
            // access neighbor
            ElementPointer neighborPointer = isIt->outside();

            int globalIdxJ = problem_.variables().indexTransport(*neighborPointer);

            // compute factor in neighbor
            GeometryType neighborGT = neighborPointer->geometry().type();
            const LocalPosition& localPosNeighbor = ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

            // neighbor cell center in global coordinates
            const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

            // take arithmetic average of absolute permeability
            permeability += soil_.K(globalPosNeighbor, *neighborPointer, localPosNeighbor);
            permeability *= 0.5;

            //get lambda_bar = lambda_n*f_w
	    if(preComput_)
        	mobBarJ=problem_.variables().mobilityWetting()[globalIdxJ]*problem_.variables().fracFlowFuncNonWetting()[globalIdxJ];
	    else
        	mobBarJ=problem_.materialLaw().mobN(1-satJ, globalPosNeighbor, *neighborPointer, localPosNeighbor)*problem_.materialLaw().fractionalW(satJ, globalPosNeighbor, *neighborPointer, localPosNeighbor);
        }
        else
        {
            //calculate lambda_n*f_w at the boundary: use a regularization with regularizationparamter eps_
            std::vector<Scalar> mobI = problem_.materialLaw().mob((eps_*satJ + (1-eps_)*satI), globalPos, element, localPos);
            std::vector<Scalar> mobJ = problem_.materialLaw().mob((eps_*satI + (1-eps_)*satJ), globalPos, element, localPos);
            mobBarI = mobI[0]*mobI[1]/(mobI[0]+mobI[1]);
            mobBarJ = mobJ[0]*mobJ[1]/(mobJ[0]+mobJ[1]);
        }

        // set result to K*grad(pc)
        FieldVector result(0);
        permeability.umv(pcGradient, result);

        // set result to f_w*lambda_n*K*grad(pc)
        result *= (mobBarI+mobBarJ)*0.5;

        return result;
    }

    CapillaryDiffusion (Problem& problem, Matrix2p<Grid, Scalar>& soil, const bool preComput = true)
    : problem_(problem), soil_(soil), preComput_(preComput), eps_(0.05)
    {}

private:
    Problem& problem_;
    Matrix2p<Grid, Scalar>& soil_;
    const bool preComput_;
    Scalar eps_;
};
}

#endif
