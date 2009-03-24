//$Id$
#ifndef DUNE_ELASTICSOIL_HH
#define DUNE_ELASTICSOIL_HH

#include <dumux/material/matrixproperties.hh>

namespace Dune
{
//! class that defines the soil parameters for the deformation of an elastic matrix
 /*! Soil parameter definition for the deformation of an elastic matrix.
  *
  *    Template parameters are:
  *
  *    - ScalarT  Floating point type used for scalars
  */
    template<class Grid, class ScalarT>
    class ElasticSoil
    {
    public:
        typedef typename Grid::Traits::template Codim<0>::Entity Element;
        typedef ScalarT Scalar;
        typedef typename Grid::ctype CoordScalar;
        enum {dim=Grid::dimension, dimWorld=Grid::dimensionworld};

        typedef Dune::FieldVector<CoordScalar,dim>      LocalPosition;
        typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

        ElasticSoil()
            {}


        const FieldVector<Scalar,2> lameParams(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const
            {
            // example for Lame parameters
                            FieldVector<Scalar,2> param;
                            param[0] = 0.443; // lambda (compressibility)
                            param[1] = 0.01; // mu (rigidity)

                            return param;
            }
     };
}
#endif
