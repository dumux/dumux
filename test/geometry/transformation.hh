//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Common
 * \brief Create a transformation function (translation + rotation)
 */
#ifndef DUMUX_TEST_GEOMETRY_TRANSFORMATION_HH
#define DUMUX_TEST_GEOMETRY_TRANSFORMATION_HH

#include <iostream>

#include <dune/common/classname.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Create a transformation function (translation) in 1D
 */
template<class ctype>
auto make1DTransformation(const ctype scale,
                          const Dune::FieldVector<ctype, 1>& translate,
                          bool verbose = true)
{
    if (verbose)
        std::cout << "Created 1D transformation with"
                  << " ctype: " << Dune::className<ctype>()
                  << ", scaling: " << scale
                  << ", translation: " << translate << std::endl;

    return [=](Dune::FieldVector<ctype, 1> p){
        p *= scale;
        p.axpy(scale, translate);
        return p;
    };
}

/*!
 * \ingroup Common
 * \brief Create a transformation function (translation + rotation) in 2D
 * \note https://en.wikipedia.org/wiki/Rotation_matrix
 */
template<class ctype>
auto make2DTransformation(const ctype scale,
                          const Dune::FieldVector<ctype, 2>& translate,
                          const ctype rotationAngle,
                          bool verbose = true)
{
    if (verbose)
        std::cout << "Created 2D transformation with"
                  << " ctype: " << Dune::className<ctype>()
                  << ", scaling: " << scale
                  << ", translation: " << translate
                  << ", rotationAngle: " << rotationAngle << std::endl;

    using std::sin; using std::cos;
    const ctype sinAngle = sin(rotationAngle);
    const ctype cosAngle = cos(rotationAngle);
    return [=](Dune::FieldVector<ctype, 2> p){
        auto tp = p;
        tp[0] = p[0]*cosAngle-p[1]*sinAngle;
        tp[1] = p[0]*sinAngle+p[1]*cosAngle;
        tp *= scale;
        return tp.axpy(scale, translate);
    };
}

/*!
 * \ingroup Common
 * \brief Create a transformation function (translation + rotation) in 3D
 * \note https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
 */
template<class ctype>
auto make3DTransformation(const ctype scale,
                          const Dune::FieldVector<ctype, 3>& translate,
                          const Dune::FieldVector<ctype, 3>& rotationAxis,
                          const ctype rotationAngle,
                          bool verbose = true)
{
    if (verbose)
        std::cout << "Created 3D transformation with"
                  << " ctype: " << Dune::className<ctype>()
                  << ", scaling: " << scale
                  << ", translation: " << translate
                  << ", rotationAxis: " << rotationAxis
                  << ", rotationAngle: " << rotationAngle << std::endl;

    using std::sin; using std::cos;
    const ctype sinAngle = sin(rotationAngle);
    const ctype cosAngle = cos(rotationAngle);
    return [=](Dune::FieldVector<ctype, 3> p){
        auto tp = p;
        tp *= cosAngle;
        tp.axpy(sinAngle, Dumux::crossProduct({rotationAxis}, p));
        tp.axpy((1.0-cosAngle)*(rotationAxis*p), rotationAxis);
        tp *= scale;
        return tp.axpy(scale, translate);
    };
}

} // end namespace Dumux

#endif
