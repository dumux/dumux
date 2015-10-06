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
    GnuplotInterface(bool persist = true) :
        plotStyle_("lines"), pipe_(0),
        fileName_(0), plotOptions_(0), plotName_(0),
        interaction_(true),
        xRangeMin_(1e100), xRangeMax_(-1e100),
        yRangeMin_(1e100), yRangeMax_(-1e100),
        xLabel_(""), yLabel_(""),
        options_(""),
        datafileSeparator_(' '),
        gnuplotPath_(GNUPLOT_EXECUTABLE),
        terminalType_("wxt")
    {
        if (persist)
            pipe_ = popen((gnuplotPath_ + " -persist").c_str(), "w"); // "w" - writing
        else
            pipe_ = popen((gnuplotPath_).c_str(), "w");
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
     * \param plottingWindowNumber Change the number of the window in which the plot is shown
     */
    void plot(const std::string &title, const unsigned int plottingWindowNumber = 0)
    {
        // set correct terminal and general options
        std::string plot = "reset\n";
        plot += "set datafile separator \'" + convertToString(datafileSeparator_) + "\'\n";
        plot += "set term wxt " + convertToString(plottingWindowNumber) + " \n";
        if (xRangeMin_ < 1e100 || xRangeMax_ > -1e100)
        {
            plot += "set xrange [" + convertToString(xRangeMin_)
                    + ":" + convertToString(xRangeMax_) + "]" + "\n";
        }
        if (yRangeMin_ < 1e100 || yRangeMax_ > -1e100)
        {
            plot += "set yrange [" + convertToString(yRangeMin_)
                    + ":" + convertToString(yRangeMax_) + "]" + "\n";
        }
        plot += "set xlabel \"" + xLabel_ + "\"\n";
        plot += "set ylabel \"" + yLabel_ + "\"\n";

        // set user defined options
        plot += options_ + "\n";

        // plot curves
        plot += "plot";
        for (unsigned int i = 0; i < fileName_.size(); ++i)
        {
            plot += + " " + fileName_[i]
                    + " " + plotOptions_[i]
                    + " title '" + plotName_[i] + "'";
            if (i < fileName_.size()-1)
                plot += ",\\";
            plot += "\n";
        }

        // live plot of the results if gnuplot is installed
#ifdef HAVE_GNUPLOT
        if (interaction_)
            executeGnuplot(plot.c_str());
#endif

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
     */
    void reset()
    {
        fileName_.resize(0);
        plotOptions_.resize(0);
        plotName_.resize(0);
        options_ = "";
    }

    /*!
     * \brief Adds a function to list of plotted lines for specific window number
     *
     * \param function Function to be plotted
     * \param plotName The name of the data set
     * \param plotOptions Specific gnuplot options passed to this plot
     */
    void addFunctionToPlot(const std::string function,
                           const std::string plotName,
                           const std::string plotOptions = "with lines")
    {
        fileName_.push_back(function);
        plotOptions_.push_back(plotOptions);
        plotName_.push_back(plotName);
    }

    /*!
     * \brief Adds a file to list of plotted lines for specific window number
     *
     * \param file Function to be plotted
     * \param plotName The name of the data set
     * \param plotOptions Specific gnuplot options passed to this plot
     */
    void addFileToPlot(const std::string file,
                       const std::string plotName,
                       const std::string plotOptions = "with lines")
    {
        fileName_.push_back("'" + file + "'");
        plotOptions_.push_back(plotOptions);
        plotName_.push_back(plotName);
    }

    /*!
     * \brief Adds a data set for specific window number and writes a data file
     *
     * \param x Vector containing the x-axis data points
     * \param y Vector containing the y-axis data points
     * \param plotName The name of the data set
     * \param plotOptions Specific gnuplot options passed to this plot
     */
    void addDataSetToPlot(const std::vector<Scalar>& x,
                          const std::vector<Scalar>& y,
                          const std::string plotName,
                          const std::string plotOptions = "with lines")
    {
        assert(x.size() == y.size());

        //write data to file
        std::string fileName = plotName + ".dat";
        std::ofstream file;
        file.open(fileName);
        for (unsigned int i = 0; i < x.size(); i++)
        {
            checkNumber(x[i], "x[i] i=" + convertToString(i) + " in " + fileName);
            checkNumber(y[i], "y[i] i=" + convertToString(i) + " in " + fileName);
            file << x[i] << datafileSeparator_ << y[i] << std::endl;
        }
        file.close();

        // adding file to list of plotted lines
        fileName_.push_back("'" + fileName + "'");
        plotOptions_.push_back(plotOptions);
        plotName_.push_back(plotName);
    }

    /*!
     * \brief Set whether a gnuplot window will be opened or not
     *
     * \param interaction Open window or not
     */
    void setInteraction(bool interaction)
    {
        interaction_ = interaction;
    }

    /*!
     * \brief Sets the label for the x-axis
     *
     * \param label The label of the x-axis
     */
    void setXlabel(const std::string& label)
    {
        xLabel_ = label;
    }

    /*!
     * \brief Sets the label for the y-axis
     *
     * \param label The label of the y-axis
     */
    void setYlabel(const std::string& label)
    {
        yLabel_ = label;
    }

    /*!
     * \brief Sets the range for the x-axis
     *
     * \param lowerEnd The lowest plotted value for the x-axis
     * \param upperEnd The highest plotted value for the x-axis
     */
    void setXRange(Scalar lowerEnd, Scalar upperEnd)
    {
        xRangeMin_ = std::min(xRangeMin_, lowerEnd);
        xRangeMax_ = std::max(xRangeMax_, upperEnd);
    }

    /*!
     * \brief Sets the range for the y-axis
     *
     * \param lowerEnd The lowest plotted value for the y-axis
     * \param upperEnd The highest plotted value for the y-axis
     */
    void setYRange(Scalar lowerEnd, Scalar upperEnd)
    {
        yRangeMin_ = std::min(yRangeMin_, lowerEnd);
        yRangeMax_ = std::max(yRangeMax_, upperEnd);
    }

    /*!
     * \brief Sets additional user-defined options
     *
     * \param option Additional line of option in gnuplot language
     */
    void setOption(std::string option)
    {
        options_ += option + "\n";
    }

    /*!
     * \brief Sets the datafile separator
     *
     * \param separator The separator sign between two data columns
     */
    void setDatafileSeparator(char separator)
    {
        datafileSeparator_ = separator;
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
    StringVector fileName_;
    StringVector plotOptions_;
    StringVector plotName_;
    bool interaction_;
    Scalar xRangeMin_;
    Scalar xRangeMax_;
    Scalar yRangeMin_;
    Scalar yRangeMax_;
    std::string xLabel_;
    std::string yLabel_;
    std::string options_;
    char datafileSeparator_;
    std::string gnuplotPath_;
    std::string terminalType_;
};
} // end of namespace
#endif // DUMUX_GNUPLOT_INTERFACE_HH
