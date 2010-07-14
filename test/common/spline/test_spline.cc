// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \brief This is a program to test the polynomial spline interpolation.
 *
 * It just prints some function to stdout. You can look at the result
 * using the following commands:
 *
----------- snip -----------
./test_spline > spline.csv
gnuplot

gnuplot> plot "spline.csv" using 1:2 w l ti "Curve", \
              "spline.csv" using 1:3 w l ti "Derivative", \
              "spline.csv" using 1:4 w p ti "Monotonical"
----------- snap -----------


*/
#include "config.h"

#include <dumux/common/spline.hh>

#include <dune/common/fvector.hh>

template <int numSamples>
void plot(bool reallyPlot)
{
    if (!reallyPlot)
        return;

    const int n = numSamples - 1;
    typedef Dune::FieldVector<double, numSamples> FV;

    double x_[] = { 0, 5, 7.5, 8.75, 9.375 };
    double y_[] = { 10, 0, 10, 0, 10 };
    double m1 =  10;
    double m2 = -10;
    FV &xs = *reinterpret_cast<FV*>(x_);
    FV &ys = *reinterpret_cast<FV*>(y_);

    Dumux::Spline<double, numSamples> sp(xs, ys, m1, m2);

    sp.printCSV(-0.1*(x_[n] - x_[0]) + x_[0],
                0.1*(x_[n] - x_[0]) + x_[n],
                1000);

    std::cerr << "Spline is monotonic: " << sp.monotonic(x_[0], x_[n]) << "\n";
};

int main(int argc, char** argv)
{
    plot<2>(false);
    plot<5>(true);
    return 0;
}
