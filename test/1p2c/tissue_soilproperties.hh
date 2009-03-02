/// $Id$
#ifndef TISSUE_SOILPROPERTIES
#define TISSUE_SOILPROPERTIES

#include <dumux/material/property_baseclasses.hh>

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class TissueSoil: public Matrix2p<Grid, Scalar>
{
public:
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::ctype CoordScalar;
    enum
        {    dim=Grid::dimension, numEq=2};

    const FieldMatrix<CoordScalar,dim,dim> &K (const FieldVector<CoordScalar,dim>& globalPos, const Element& element, const FieldVector<CoordScalar,dim>& localPos)
    {
        if (10<globalPos[0] && globalPos[0]<12 && 10<globalPos[1] && globalPos[1]<12)
            return permloc_;                        //tumor tissue
        else
            return permlocWell_;                    //healthy tissue
    }

    double porosity(const FieldVector<CoordScalar,dim>& globalPos, const Element& element, const FieldVector<CoordScalar,dim>& localPos) const
    {
        if (10<globalPos[0] && globalPos[0]<12 && 10<globalPos[1] && globalPos[1]<12)
            return 0.31;                        //tumor tissue
        else
            return 0.13;                    //healthy tissue
    }

    double tortuosity (const FieldVector<CoordScalar,dim>& globalPos, const Element& element, const FieldVector<CoordScalar,dim>& localPos) const
    {
        if (10<globalPos[0] && globalPos[0]<12 && 10<globalPos[1] && globalPos[1]<12)
            return 0.706;                        //tumor tissue
        else
            return 0.280;                    //healthy tissue
    }

    virtual double Sr_w(const FieldVector<CoordScalar,dim>& globalPos, const Element& element, const FieldVector<CoordScalar,dim>& localPos, const double T = 283.15) const
    {
        return 0;
    }

    virtual double Sr_n(const FieldVector<CoordScalar,dim>& globalPos, const Element& element, const FieldVector<CoordScalar,dim>& localPos, const double T = 283.15) const
    {
        return 0;
    }

    virtual std::vector<double> paramRelPerm(const FieldVector<CoordScalar,dim>& globalPos, const Element& element, const FieldVector<CoordScalar,dim>& localPos, const double T = 283.15) const
    {
        std::vector<double> param(0);
        return param;
    }

    TissueSoil()
        :Matrix2p<Grid,Scalar>()
    {
        permloc_ = 0;
        permlocWell_ = 0;
        for (int i = 0; i < dim; i++)
            permloc_[i][i] = 2.142e-11;             //tumor tissue
        for (int i = 0; i < dim; i++)
            permlocWell_[i][i] = 4.424e-12;        //healthy tissue
    }

private:
    Dune::FieldMatrix<CoordScalar,dim,dim> permloc_, permlocWell_;
};
} // end namespace
#endif
