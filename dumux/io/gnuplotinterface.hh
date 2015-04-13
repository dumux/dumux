// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Interface for passing data sets to a file and plotting them, if gnuplot
 *        is installed.
 *
 * The data sets for a specific window have to be passed by the addDataSet function
 * and then plotted by using the plot function.
 */

#ifndef HAVE_GNUPLOT
#warning Gnuplot has not been found by CMake, no output possible.
#define GNUPLOT_EXECUTABLE "/usr/bin/gnuplot"
#endif

#ifndef DUMUX_GNUPLOT_INTERFACE_HH
#define DUMUX_GNUPLOT_INTERFACE_HH

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/stdstreams.hh>

namespace Dumux
{
/*!
 * \brief Interface for passing data sets to a file and plotting them, if gnuplot
 *        is installed.
 */
template<class Scalar>
class GnuplotInterface
{
public:
    typedef std::vector<std::string> StringVector;
    enum PlotStyle
    {
        lines, points, linesPoints, impulses, dots
    };

    GnuplotInterface() :
        plotStyle_("lines"), pipe_(0),
        fileName_(0), plotName_(0),
        xRangeMin_(0), xRangeMax_(0),
        yRangeMin_(0), yRangeMax_(0),
        xLabel_(0), yLabel_(0),
        gnuplotPath_(GNUPLOT_EXECUTABLE),
        terminalType_("wxt"),
        numberOfGnuplotWindows_(8+1)
    {
        {
            pipe_ = popen((gnuplotPath_ + " -persist").c_str(), "w"); // "w" - writing
        }
        fileName_.resize(numberOfGnuplotWindows_);
        plotName_.resize(numberOfGnuplotWindows_);
        xRangeMin_.resize(numberOfGnuplotWindows_, 1e100);
        xRangeMax_.resize(numberOfGnuplotWindows_, -1e100);
        yRangeMin_.resize(numberOfGnuplotWindows_, 1e100);
        yRangeMax_.resize(numberOfGnuplotWindows_, -1e100);
        xLabel_.resize(numberOfGnuplotWindows_, "");
        yLabel_.resize(numberOfGnuplotWindows_, "");
    }

    ~GnuplotInterface()
    {
        if (pclose(pipe_) == -1)
            assert("Could not close pipe to Gnuplot!");
    }
    
    /*!
     * \brief Plots the files for a specific window number, writes a gnuplot and png file.
     *
     * \param title The name of the output file
     * \param plottingWindowNumber The ID of the specific plot
     */
    void plot(const std::string &title,
              const unsigned int plottingWindowNumber,
              bool interaction)
    {
        // setting correct terminal
        std::string plot = "set term ";

        if (interaction)
            plot += "wxt " + numberToString(plottingWindowNumber) + "\n";
        else
            plot += " pngcairo size 800,600 \nset output \"" + title + ".png\"\n";

        plot += "set yrange [" + numberToString(yRangeMin_[plottingWindowNumber])
                + ":" + numberToString(yRangeMax_[plottingWindowNumber]) + "]" + "\n";
        plot += "set xrange [" + numberToString(xRangeMin_[plottingWindowNumber])
                + ":" + numberToString(xRangeMax_[plottingWindowNumber]) + "]" + "\n";

        plot += "set xlabel \"" + xLabel_[plottingWindowNumber] + "\"\n";
        plot += "set ylabel \"" + yLabel_[plottingWindowNumber] + "\"\n";
        
        // plot curves
        plot += "plot";
        for (unsigned int i = 0; i < fileName_[plottingWindowNumber].size(); ++i)
        {
            plot += + " " + fileName_[plottingWindowNumber][i]
                    + " title " + plotName_[plottingWindowNumber][i]
                    + " with " + plotStyle_
                    + ",";
        }
        plot[plot.size()-1] = '\n'; // remove last ","

        if (!interaction)
        {
            std::string fileName = title + ".gp";
            std::ofstream file;
            file.open(fileName);
            file << plot;
            file.close();
        }

#ifdef HAVE_GNUPLOT
        executeGnuplot(plot.c_str());
#endif
    }

    /*!
     * \brief Deletes all plots from a plotting window
     *
     * \param plottingWindowNumber The ID of the specific plot
     */
    void reset(const unsigned int plottingWindowNumber)
    {
        fileName_[plottingWindowNumber].resize(0);
        plotName_[plottingWindowNumber].resize(0);
    }

    /*!
     * \brief Adds a function to list of plotted lines for specific window number
     *
     * \param function Function to be plotted
     * \param plotName The name of the data set
     * \param plottingWindowNumber The ID of the specific plot
     */
    void addFunctionToPlot(const std::string &function,
                           const std::string &plotName,
                           const unsigned int plottingWindowNumber)
    {
        fileName_[plottingWindowNumber].push_back(function);
        plotName_[plottingWindowNumber].push_back("'" + plotName + "'");
    }

    /*!
     * \brief Adds a file to list of plotted lines for specific window number
     *
     * \param addFile Function to be plotted
     * \param plotName The name of the data set
     * \param plottingWindowNumber The ID of the specific plot
     */
    void addFileToPlot(const std::string &file,
                       const std::string &plotName,
                       const unsigned int plottingWindowNumber)
    {
        fileName_[plottingWindowNumber].push_back("'" + file + "'");
        plotName_[plottingWindowNumber].push_back("'" + plotName + "'");
    }

    /*!
     * \brief Adds a data set for specific window number and writes a data file
     *
     * \param x Vector containing the x-axis data points
     * \param y Vector containing the y-axis data points
     * \param plotName The name of the data set
     * \param plottingWindowNumber The ID of the specific plot
     */
    void addDataSetToPlot(const std::vector<Scalar>& x,
                          const std::vector<Scalar>& y,
                          const std::string &plotName,
                          const unsigned int plottingWindowNumber)
    {
        assert(x.size() == y.size());

        //write data to file
        std::string fileName = numberToString(plottingWindowNumber) + "_" + plotName + ".dat";
        std::ofstream file;
        file.open(fileName);
        for (unsigned int i = 0; i < x.size(); i++)
        {
            checkNumber(x[i], "x[i] i=" + numberToString(i) + " in " + fileName);
            checkNumber(y[i], "y[i] i=" + numberToString(i) + " in " + fileName);
            file << x[i] << " " << y[i] << std::endl;
        }
        file.close();

        // adding file to list of plotted lines
        fileName_[plottingWindowNumber].push_back("'" + fileName + "'");
        plotName_[plottingWindowNumber].push_back("'" + fileName + "'");
    }

    /*!
     * \brief Sets the label for the x-axis
     *
     * \param label The label of the x-axis
     * \param plottingWindowNumber The ID of the specific plot
     */
    void setXlabel(const std::string& label,
                   const unsigned int plottingWindowNumber)
    {
        xLabel_[plottingWindowNumber] = label;
    }

    /*!
     * \brief Sets the label for the y-axis
     *
     * \param label The label of the y-axis
     * \param plottingWindowNumber The ID of the specific plot
     */
    void setYlabel(const std::string& label,
                   const unsigned int plottingWindowNumber)
    {
        yLabel_[plottingWindowNumber] = label;
    }

    /*!
     * \brief Sets the range for the x-axis
     *
     * \param lowerEnd The lowest plotted value for the x-axis
     * \param upperEnd The highest plotted value for the x-axis
     * \param plottingWindowNumber The ID of the specific plot
     */
    void setXRange(Scalar lowerEnd,
                   Scalar upperEnd,
                   const unsigned int plottingWindowNumber)
    {
        xRangeMin_[plottingWindowNumber] = std::min(xRangeMin_[plottingWindowNumber], lowerEnd);
        xRangeMax_[plottingWindowNumber] = std::max(xRangeMax_[plottingWindowNumber], upperEnd);
    }

    /*!
     * \brief Sets the range for the y-axis
     *
     * \param lowerEnd The lowest plotted value for the y-axis
     * \param upperEnd The highest plotted value for the y-axis
     * \param plottingWindowNumber The ID of the specific plot
     */
    void setYRange(Scalar lowerEnd,
                   Scalar upperEnd,
                   const unsigned int plottingWindowNumber)
    {
        yRangeMin_[plottingWindowNumber] = std::min(yRangeMin_[plottingWindowNumber], lowerEnd);
        yRangeMax_[plottingWindowNumber] = std::max(yRangeMax_[plottingWindowNumber], upperEnd);
    }

    /*!
     * \brief Sets the plotting style for the data sets
     *
     * \param style Plot style of the data sets
     */
    void setStyle(const PlotStyle& style)
    {
        switch (style)
        {
        case lines:
            plotStyle_ = "lines";
            break;
        case points:
            plotStyle_ = "points";
            break;
        case linesPoints:
            plotStyle_ = "linespoints";
            break;
        case impulses:
            plotStyle_ = "impulses";
            break;
        case dots:
            plotStyle_ = "dots";
            break;
        default:
            assert(!"Unknown plot style");
        }
    }

private:
    // Give plot command to gnuplot
    void executeGnuplot(const std::string& plotCommand) const
    {
        fputs((plotCommand + "\n").c_str(), pipe_);
        fflush(pipe_);
    }

    // Check validity of number
    void checkNumber(Scalar number, std::string text = "") const
    {
        if (std::isnan(number))
            Dune::dwarn << "warning: " << text << " is not a number, adjust your data range" << std::endl;
        if (std::isinf(number))
            Dune::dwarn << "warning: " << text << " is infinity, adjust your data range" << std::endl;
    }

    // Convert number to string
    template<class T> std::string numberToString(const T number) const
    {
        std::ostringstream stream;
        stream << number;
        return stream.str();
    }

    std::string plotStyle_;
    std::FILE * pipe_;
    std::vector<StringVector> fileName_;
    std::vector<StringVector> plotName_;
    std::vector<Scalar> xRangeMin_;
    std::vector<Scalar> xRangeMax_;
    std::vector<Scalar> yRangeMin_;
    std::vector<Scalar> yRangeMax_;
    StringVector xLabel_;
    StringVector yLabel_;
    std::string gnuplotPath_;
    std::string terminalType_;
    const unsigned int numberOfGnuplotWindows_;
};
} // end of namespace
#endif // DUMUX_GNUPLOT_INTERFACE_HH
