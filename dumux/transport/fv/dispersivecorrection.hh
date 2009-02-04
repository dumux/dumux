// $Id$
#ifndef DUNE_DISPERSIVECORRECTION_HH
#define DUNE_DISPERSIVECORRECTION_HH

#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/fractionalflow/fractionalflowproblem.hh"
//#include "dumux/material/property_baseclasses.hh"

namespace Dune
{

template<class Grid, class Scalar, class VC>
class DispersiveCorrection: public DiffusivePart<Grid, Scalar>
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
    virtual Dune::FieldVector<Scalar, dim> operator() (const Element& element, const int numberInSelf,
            const Scalar satIntersection, const Dune::FieldVector<Scalar, dim>& satGradient, const Scalar time,
            const Scalar satI, const Scalar satJ) const
    {
        Scalar result = 0;
        Scalar helpresult = 0;

        int globalIdxI = problem_.variables.transMapper.map(element);

        int colNumI = problem_.soil.getDispersion()[globalIdxI].size();

        IntersectionIterator endis = element.ilevelend();
        for (IntersectionIterator isIt = element.ilevelbegin(); isIt != endis; ++isIt)
        {
            if(isIt->numberInSelf() != numberInSelf)
            {
                continue;
            }

            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();

                int globalIdxJ = problem_.variables.transMapper.map(*neighborPointer);

                int colNumJ = problem_.soil.getDispersion()[globalIdxJ].size();

                for (int i=0;i<colNumI;i++)
                {
                    if (problem_.soil.getDispersionSat()[globalIdxI][i]> satI)
                    {
                        Scalar satdiff1 = problem_.soil.getDispersionSat()[globalIdxI][i] - satI;
                        Scalar satdiff2 = 1e100;
                        if (i)
                        {
                            satdiff2 = satI - problem_.soil.getDispersionSat()[globalIdxI][i-1];
                        }
                        if (satdiff1 < satdiff2)
                        {
                            result = problem_.soil.getDispersion()[globalIdxI][i]*satGradient;
//                                                    std::cout<<"result1 = "<<problem_.soil.getDispersion()[globalIdxI][i]<<std::endl;
                        }
                        else
                        {
                            result = problem_.soil.getDispersion()[globalIdxI][i-1]*satGradient;
//                                                    std::cout<<"result2 = "<<problem_.soil.getDispersion()[globalIdxI][i-1]<<std::endl;
                            //                        std::cout<<"satGrad = "<<satGradient<<std::endl;
                        }
                        break;
                    }
                    if (i==(colNumI-1))
                    {
                        result = problem_.soil.getDispersion()[globalIdxI][i]*satGradient;
                        break;
                    }
                }
                for (int i=0;i<colNumJ;i++)
                {
                    if (problem_.soil.getDispersionSat()[globalIdxJ][i]> satJ)
                    {
                        Scalar satdiff1 = problem_.soil.getDispersionSat()[globalIdxJ][i] - satJ;
                        Scalar satdiff2 = 1e100;
                        if (i)
                        {
                            satdiff2 = satJ - problem_.soil.getDispersionSat()[globalIdxJ][i-1];
                        }
                        if (satdiff1 < satdiff2)
                        {
                            helpresult = problem_.soil.getDispersion()[globalIdxJ][i]*satGradient;
//                                                    std::cout<<"helpresult1 = "<<problem_.soil.getDispersion()[globalIdxJ][i]<<std::endl;
                        }
                        else
                        {
                            //                        std::cout<<"D = "<<problem_.soil.getDispersion()[globalIdxJ][i-1]<<std::endl;
                            //                        std::cout<<"satG = "<<satGradient<<std::endl;
                            helpresult = problem_.soil.getDispersion()[globalIdxJ][i-1]*satGradient;
//                                                    std::cout<<"helpresult2 = "<<problem_.soil.getDispersion()[globalIdxJ][i-1]<<std::endl;
                        }
                        break;
                    }
                    if (i==(colNumJ-1))
                    {
                        helpresult = problem_.soil.getDispersion()[globalIdxJ][i]*satGradient;
                        break;
                    }
                }

                if (result*helpresult == 0)
                {
                    result=0;
                }
                else
                {
                    result = 0.5*(result+helpresult);
//                    std::cout<<"result 1 = "<<result<<std::endl;
                }
            }
            else
            {
                for (int i=0;i<colNumI;i++)
                {
                    if (problem_.soil.getDispersionSat()[globalIdxI][i]> satJ)
                    {
                        Scalar satdiff1 = problem_.soil.getDispersionSat()[globalIdxI][i] - satJ;
                        Scalar satdiff2 = 1e100;
                        if (i)
                        {
                            satdiff2 = satJ - problem_.soil.getDispersionSat()[globalIdxI][i-1];
                        }
                        if (satdiff1 < satdiff2)
                        {
                            result = problem_.soil.getDispersion()[globalIdxI][i]*satGradient;
//                                                    std::cout<<"result1Bound = "<<problem_.soil.getDispersion()[globalIdxI][i]<<std::endl;
                        }
                        else
                        {
                            result = problem_.soil.getDispersion()[globalIdxI][i-1]*satGradient;
//                                                    std::cout<<"result2Bound = "<<problem_.soil.getDispersion()[globalIdxI][i]<<std::endl;
                        }
                        break;
                    }
                    if (i==(colNumI-1))
                    {
                        result = problem_.soil.getDispersion()[globalIdxI][i]*satGradient;
                        break;
                    }
                }
            }
        }

        FieldVector<Scalar, dim> returnVector(result);
//        std::cout<<"result = "<<result<<std::endl;
//        std::cout<<"returnVector = "<<returnVector<<std::endl;
        return returnVector;

    }

    DispersiveCorrection (FractionalFlowProblem<Grid, Scalar, VC>& prob)
    : problem_(prob)

    {}

private:
    FractionalFlowProblem<Grid, Scalar, VC>& problem_;
};
}

#endif
