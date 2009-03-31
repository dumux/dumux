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
 * @author Bernd Flemisch, last changed by Markus Wolff
 */
namespace Dune
{
/*!\ingroup diffPart
 * @brief  Base class for defining the diffusive part of an advection-diffusion equation
 */
template<class Grid, class Scalar, class VC, class Problem = TransportProblem<Grid, Scalar, VC> >
class CapillaryDiffusion : public DiffusivePart<Grid,Scalar>
{
    enum{dim = Grid::dimension,dimWorld = Grid::dimensionworld};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    virtual FieldVector operator() (const Element& element, const int faceNumber,
                                    const Scalar satIntersection, const FieldVector& satGradient, const Scalar time,
                                    const Scalar satI, const Scalar satJ) const
    {
        // cell geometry type
        GeometryType gt = element.geometry().type();

        // cell center in reference element
        const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0,0);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = element.geometry().global(localPos);

        // get absolute permeability of cell
        FieldMatrix K(soil_.K(globalPos,element,localPos));

        IntersectionIterator isItEnd = element.ilevelend();
        IntersectionIterator isIt = element.ilevelbegin();
        for (; isIt != isItEnd; ++isIt)
        {
            if(isIt->numberInInside() == faceNumber)
                break;
        }

        // get geometry type of face
        GeometryType faceGT = isIt->geometry().type();

        // center in face's reference element
        const Dune::FieldVector<Scalar,dim-1>& faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

        FieldVector unitOuterNormal = isIt->unitOuterNormal(faceLocal);
        //std::cout<<"unitOuterNormaldiff"<<unitOuterNormal<<std::endl;

        //get capillary pressure gradient
        Scalar dPdSI=problem_.materialLaw().dPdS(satI,globalPos,element,localPos);
        //get lambda_bar = lambda_n*f_w
        Scalar mobBarI=problem_.materialLaw().mobN(1-satI,globalPos,element,localPos)*problem_.materialLaw().fractionalW(satI,globalPos,element,localPos);

        Scalar dPdSJ;
        Scalar mobBarJ;

        if (isIt->neighbor()) {
            // access neighbor
            ElementPointer neighborPointer = isIt->outside();

            // compute factor in neighbor
            GeometryType neighborGT = neighborPointer->geometry().type();
            const LocalPosition& localPosNeighbor = ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

            // neighbor cell center in global coordinates
            const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

            // take arithmetic average of absolute permeability
            K += soil_.K(globalPosNeighbor, *neighborPointer, localPosNeighbor);
            K *= 0.5;

            //get capillary pressure gradient
            dPdSJ=problem_.materialLaw().dPdS(satJ, globalPosNeighbor, *neighborPointer, localPosNeighbor);
            //get lambda_bar = lambda_n*f_w
            mobBarJ=problem_.materialLaw().mobN(1-satJ, globalPosNeighbor, *neighborPointer, localPosNeighbor)*problem_.materialLaw().fractionalW(satJ, globalPosNeighbor, *neighborPointer, localPosNeighbor);
        }
        else
        {
            dPdSJ = dPdSI;
            mobBarJ = mobBarI;
        }

        // set result to grad(S)
        FieldVector helpResult(satGradient);

        // set result to (dp_c/dS)*grad(S)
        helpResult *= (dPdSI+dPdSJ)*0.5;

        // set result to K*((dp_c/dS)*grad(S))
        FieldVector result(0);
        K.umv(helpResult, result);

        // set result to f_w*lambda_n*K*((dp_c/dS)*grad(S))
        result *= (mobBarI+mobBarJ)*0.5;

        return result;
    }

    CapillaryDiffusion (Problem& problem, Matrix2p<Grid, Scalar>& soil)
        : problem_(problem), soil_(soil)
    { }

private:
    Problem& problem_;
    Matrix2p<Grid, Scalar>& soil_;

};
}

#endif
