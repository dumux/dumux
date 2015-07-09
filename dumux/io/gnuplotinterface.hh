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
 *
 * \todo the use of the plottingWindowNumber in each function could be replaced, but
 *       a new GnuplotInterface object has to be created for each plot
 */

#ifndef DUMUX_GNUPLOT_INTERFACE_HH
#define DUMUX_GNUPLOT_INTERFACE_HH

#if !HAVE_GNUPLOT
#warning Gnuplot has not been found by CMake, no output possible.
#define GNUPLOT_EXECUTABLE "/usr/bin/gnuplot"
#endif

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <dune/common/deprecated.hh>
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

    //! \brief The constructor
    GnuplotInterface() :
        plotStyle_("lines"), pipe_(0),
        fileName_(0), plotOptions_(0), plotName_(0),
        xRangeMin_(0), xRangeMax_(0),
        yRangeMin_(0), yRangeMax_(0),
        xLabel_(0), yLabel_(0),
        options_(0),
        datafileSeparator_(0),
        gnuplotPath_(GNUPLOT_EXECUTABLE),
        terminalType_("wxt"),
        numberOfGnuplotWindows_(8+1)
    {
        {
            pipe_ = popen((gnuplotPath_ + " -persist").c_str(), "w"); // "w" - writing
        }
        fileName_.resize(numberOfGnuplotWindows_);
        plotOptions_.resize(numberOfGnuplotWindows_);
        plotName_.resize(numberOfGnuplotWindows_);
        xRangeMin_.resize(numberOfGnuplotWindows_, 1e100);
        xRangeMax_.resize(numberOfGnuplotWindows_, -1e100);
        yRangeMin_.resize(numberOfGnuplotWindows_, 1e100);
        yRangeMax_.resize(numberOfGnuplotWindows_, -1e100);
        xLabel_.resize(numberOfGnuplotWindows_, "");
        yLabel_.resize(numberOfGnuplotWindows_, "");
        options_.resize(numberOfGnuplotWindows_, "");
        datafileSeparator_.resize(numberOfGnuplotWindows_, ' ');
    }

    //! \brief The destructor
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
     * \param interaction Specifies whether a live output via a gnuplot window is wanted
     */
    void plot(const std::string &title,
              const unsigned int plottingWindowNumber,
              bool interaction)
    {
        // set correct terminal and general options
        std::string plot = "reset\n";
        plot += "set datafile separator \'" + convertToString(datafileSeparator_[plottingWindowNumber]) + "\'\n";
        plot += "set term wxt " + convertToString(plottingWindowNumber) + "\n";
        if (xRangeMin_[plottingWindowNumber] < 1e100 || xRangeMax_[plottingWindowNumber] > -1e100)
        {
            plot += "set xrange [" + convertToString(xRangeMin_[plottingWindowNumber])
                    + ":" + convertToString(xRangeMax_[plottingWindowNumber]) + "]" + "\n";
        }
        if (yRangeMin_[plottingWindowNumber] < 1e100 || yRangeMax_[plottingWindowNumber] > -1e100)
        {
            plot += "set yrange [" + convertToString(yRangeMin_[plottingWindowNumber])
                    + ":" + convertToString(yRangeMax_[plottingWindowNumber]) + "]" + "\n";
        }
        plot += "set xlabel \"" + xLabel_[plottingWindowNumber] + "\"\n";
        plot += "set ylabel \"" + yLabel_[plottingWindowNumber] + "\"\n";

        // set user defined options
        plot += options_[plottingWindowNumber] + "\n";

        // plot curves
        plot += "plot";
        for (unsigned int i = 0; i < fileName_[plottingWindowNumber].size(); ++i)
        {
            plot += + " " + fileName_[plottingWindowNumber][i]
                    + " " + plotOptions_[plottingWindowNumber][i]
                    + " title '" + plotName_[plottingWindowNumber][i] + "'";
            if (i < fileName_[plottingWindowNumber].size()-1)
                plot += ",\\";
            plot += "\n";
        }

        // live plot of the results if gnuplot is installed
        if (interaction)
        {
#ifdef HAVE_GNUPLOT
            executeGnuplot(plot.c_str());
#endif
        }

        // create a gnuplot file for output a png
        plot += "\n";
        plot += "set term pngcairo size 800,600\n";
        plot += "set output \"" + title + ".png\"\n";
        plot += "replot\n";
        std::string fileName = title + ".gp";
        std::ofstream file;
        file.open(fileName);
        file << plot;
        file.close();
    }

    /*!
     * \brief Deletes all plots from a plotting window and resets user-defined options
     *
     * \param plottingWindowNumber The ID of the specific plot
     */
    void reset(const unsigned int plottingWindowNumber)
    {
        fileName_[plottingWindowNumber].resize(0);
        plotOptions_[plottingWindowNumber].resize(0);
        plotName_[plottingWindowNumber].resize(0);
        options_[plottingWindowNumber].resize(0);
    }

    /*!
     * \brief Adds a function to list of plotted lines for specific window number
     *
     * \param function Function to be plotted
     * \param plotName The name of the data set
     * \param plottingWindowNumber The ID of the specific plot
     * \param plotOptions Specific gnuplot options passed to this plot
     */
    void addFunctionToPlot(const std::string function,
                           const std::string plotName,
                           const unsigned int plottingWindowNumber,
                           const std::string plotOptions = "with lines")
    {
        fileName_[plottingWindowNumber].push_back(function);
        plotOptions_[plottingWindowNumber].push_back(plotOptions);
        plotName_[plottingWindowNumber].push_back(plotName);
    }

    /*!
     * \brief Adds a file to list of plotted lines for specific window number
     *
     * \param file Function to be plotted
     * \param plotName The name of the data set
     * \param plottingWindowNumber The ID of the specific plot
     * \param plotOptions Specific gnuplot options passed to this plot
     */
    void addFileToPlot(const std::string file,
                       const std::string plotName,
                       const unsigned int plottingWindowNumber,
                       const std::string plotOptions = "with lines")
    {
        fileName_[plottingWindowNumber].push_back("'" + file + "'");
        plotOptions_[plottingWindowNumber].push_back(plotOptions);
        plotName_[plottingWindowNumber].push_back(plotName);
    }

    /*!
     * \brief Adds a data set for specific window number and writes a data file
     *
     * \param x Vector containing the x-axis data points
     * \param y Vector containing the y-axis data points
     * \param plotName The name of the data set
     * \param plottingWindowNumber The ID of the specific plot
     * \param plotOptions Specific gnuplot options passed to this plot
     */
    void addDataSetToPlot(const std::vector<Scalar>& x,
                          const std::vector<Scalar>& y,
                          const std::string plotName,
                          const unsigned int plottingWindowNumber,
                          const std::string plotOptions = "with lines")
    {
        assert(x.size() == y.size());

        //write data to file
        std::string fileName = convertToString(plottingWindowNumber) + "_" + plotName + ".dat";
        std::ofstream file;
        file.open(fileName);
        for (unsigned int i = 0; i < x.size(); i++)
        {
            checkNumber(x[i], "x[i] i=" + convertToString(i) + " in " + fileName);
            checkNumber(y[i], "y[i] i=" + convertToString(i) + " in " + fileName);
            file << x[i] << datafileSeparator_[plottingWindowNumber] << y[i] << std::endl;
        }
        file.close();

        // adding file to list of plotted lines
        fileName_[plottingWindowNumber].push_back("'" + fileName + "'");
        plotOptions_[plottingWindowNumber].push_back(plotOptions);
        plotName_[plottingWindowNumber].push_back(plotName);
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
     * \brief Sets additional user-defined options
     *
     * \param option Additional line of option in gnuplot language
     * \param plottingWindowNumber The ID of the specific plot
     */
    void setOption(std::string option,
                   const unsigned int plottingWindowNumber)
    {
        options_[plottingWindowNumber] += option + "\n";
    }

    /*!
     * \brief Sets the datafile separator
     *
     * \param separator The separator sign between two data columns
     * \param plottingWindowNumber The ID of the specific plot
     */
    void setDatafileSeparator(char separator,
                              const unsigned int plottingWindowNumber)
    {
        datafileSeparator_[plottingWindowNumber] = separator;
    }

    /*!
     * \brief Sets the plotting style for the data sets
     *
     * \param style Plot style of the data sets
     */
    DUNE_DEPRECATED_MSG("setStyle() functionality has been replaced by directly passing gnuplot options to the add..ToPlot function")
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

    // Convert number or character to string
    template<class T>
    std::string convertToString(const T input) const
    {
        std::ostringstream stream;
        stream << input;
        return stream.str();
    }

    std::string plotStyle_;
    std::FILE * pipe_;
    std::vector<StringVector> fileName_;
    std::vector<StringVector> plotOptions_;
    std::vector<StringVector> plotName_;
    std::vector<Scalar> xRangeMin_;
    std::vector<Scalar> xRangeMax_;
    std::vector<Scalar> yRangeMin_;
    std::vector<Scalar> yRangeMax_;
    StringVector xLabel_;
    StringVector yLabel_;
    StringVector options_;
    std::vector<char> datafileSeparator_;
    std::string gnuplotPath_;
    std::string terminalType_;
    const unsigned int numberOfGnuplotWindows_;
};
} // end of namespace
#endif // DUMUX_GNUPLOT_INTERFACE_HH
