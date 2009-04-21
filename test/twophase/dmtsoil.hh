#ifndef DUNE_DMTSOIL_HH
#define DUNE_DMTSOIL_HH

/**
 * @file
 * @brief  Class for defining an instance of a Matrix2p soil
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class DMTSoil: public Matrix2p<Grid,Scalar>
{
public:
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename Grid::ctype DT;
    enum {dim=Grid::dimension, numEq=1};

    // define PERMEABILITY tensor
    virtual const FieldMatrix<DT,dim,dim> &K (const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi) const
    {
        FieldMatrix<DT,dim,dim> result(0);

        std::vector<double>& parameters = gridPtr_.parameters(element);
        switch ((int)parameters[0])
        {
        case 1:
            result[0][0] = result[1][1] = 1e-10;
            break;
        case 2:
            result[0][0] = result[1][1] = 5e-12;
            break;
        case 3:
            result[0][0] = result[1][1] = 5e-12;
            break;
        case 4:
            if (x[1] >= 0.55)
            {
                if ( x[0] >= 0.2)
                    result[0][0] = result[1][1] = 1.0e-19;
                else
                    result[0][0] = result[1][1] = 1.0e-7;
            }
            else
                result[0][0] = result[1][1] = 5.0e-12;
            break;
        }

        return result;
    }

    virtual double porosity(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi) const
    {
        std::vector<double>& parameters = gridPtr_.parameters(element);

        switch ((int)parameters[0])
        {
        case 1:
        case 2:
        case 3:
            return 0.3;
            break;
        case 4:
            return 0.25;
            break;
        }

        return 0.0;
    }

    virtual double Sr_w(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T) const
    {
        return 0.0;
    }

    virtual double Sr_n(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T) const
    {
        return 0.1;
    }

    virtual typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi) const
    {
        return Matrix2p<Grid,Scalar>::linear;
    }

    virtual std::vector<double> paramRelPerm(const FieldVector<DT,dim>& x, const Entity& element, const FieldVector<DT,dim>& xi, const double T) const
    {
        std::vector<double> param(4);

        param[0] = param[1] = 0.0; // parameters for capillary pressure (= 0)
        param[2] = 0.0; // constant k_rw
        param[3] = 0.9; // constant k_rn

        return param;
    }

    DMTSoil(GridPtr<Grid>& gP)
        : Matrix2p<Grid,Scalar>(), gridPtr_(gP)
    {
    }

    ~DMTSoil()
    {}

private:
    GridPtr<Grid>& gridPtr_;
};

} // end namespace
#endif

