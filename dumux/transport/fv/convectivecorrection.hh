// $Id$
#ifndef DUNE_CONVECTIVECORRECTION_HH
#define DUNE_CONVECTIVECORRECTION_HH

#include "dumux/transport/fv/convectivepart.hh"
#include "dumux/fractionalflow/fractionalflowproblem.hh"

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar, class VC>
class ConvectiveCorrection: public ConvectivePart<Grid, Scalar>
{
    enum
        {
            dim = Grid::dimension, dimWorld = Grid::dimensionworld
        };
    typedef    typename Grid::LevelGridView GridView;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    virtual double operator() (const Element& element, const Scalar sat, const GlobalPosition& globalPosFace) const
    {
        double result = 0;

        int globalIdx = problem.variables.transMapper.map(element);

        int colNum = problem.soil.getM()[globalIdx].size();

        // element geometry
        const Geometry& geometry = element.geometry();

        GeometryType gt = element.type();
        // cell center in reference element
        const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = geometry.global(localPos);

        IntersectionIterator isItEnd = element.ilevelend();
        for (IntersectionIterator isIt = element.ilevelbegin(); isIt != isItEnd; ++isIt)
        {
            // get geometry type of face
            GeometryType faceGT = isIt->intersectionSelfLocal().type();

            // center in face's reference element
            const FieldVector<Scalar,dim-1>& faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

            // center of face in global coordinates
            const GlobalPosition& globalPosFaceCheck = isIt->intersectionGlobal().global(faceLocal);

            if (globalPosFace[0] <= globalPosFaceCheck[0] + eps_ && globalPosFace[0] >= globalPosFaceCheck[0] - eps_
                && globalPosFace[1] <= globalPosFaceCheck[1] + eps_ && globalPosFace[1] <= globalPosFaceCheck[1] + eps_)
            {
                int faceNumber = isIt->numberInSelf();
                for (int i=0;i<colNum;i++)
                {
                    if (sat <= problem.soil.Sr_w(globalPos, element, localPos))
                    {
                        return 0;
                    }
                    if (problem.soil.getMSat()[globalIdx][i][faceNumber]>= sat)
                    {
                        //                        std::cout<<"sat = "<<sat<<"satM = "<<problem.soil.getMSat()[globalIdx][i][faceNumber]<<std::endl;
                        double satdiff1 = problem.soil.getMSat()[globalIdx][i][faceNumber] - sat;
                        double satdiff2 = 1e100;
                        if (i)
                        {
                            satdiff2 = sat - problem.soil.getMSat()[globalIdx][i-1][faceNumber];
                        }
                        if (satdiff1 < satdiff2)
                        {
                            result = problem.soil.getM()[globalIdx][i][faceNumber];
                            //                            std::cout<<"m = "<<problem.soil.getM()[globalIdx][i][faceNumber]<<std::endl;
                            //                            std::cout<<"sat = "<<sat<<"satM = "<<problem.soil.getMSat()[globalIdx][i][faceNumber]<<std::endl;
                        }
                        else
                        {
                            result = problem.soil.getM()[globalIdx][i-1][faceNumber];
                            //                            std::cout<<"m = "<<problem.soil.getM()[globalIdx][i][faceNumber]<<std::endl;
                            //                            std::cout<<"sat = "<<sat<<"satM = "<<problem.soil.getMSat()[globalIdx][i][faceNumber]<<std::endl;

                        }
                        break;
                    }
                    if (i==(colNum-1))
                    {
                        result = problem.soil.getM()[globalIdx][i][faceNumber];
                        break;
                    }
                }
            }
        }

        return result;

    }

    ConvectiveCorrection (FractionalFlowProblem<Grid, Scalar, VC>& prob)
        : problem(prob), eps_(1e-6)

    {}

private:
    FractionalFlowProblem<Grid, Scalar, VC>& problem;
    Scalar eps_;
};
}

#endif
