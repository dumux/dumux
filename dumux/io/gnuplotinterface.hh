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
    GnuplotInterface(bool persist = true) :
        pipe_(0), openPlotWindow_(true), persist_(persist),
        terminalType_("x11"),
        datafileSeparator_(' '), linetype_("solid"),
        xRangeIsSet_(false), yRangeIsSet_(false),
        xLabel_(""), yLabel_(""),
        gnuplotPath_(GNUPLOT_EXECUTABLE)
    {
        open(persist_);
        resetPlot();
    }

    //! \brief The destructor
    ~GnuplotInterface()
    {
        close();
    }

    /*!
     * \brief Plots the files for a specific window number, writes a gnuplot and png file.
     *
     * \param title The name of the output file
     * \param plottingWindowNumber Change the number of the window in which the plot is shown
     * \param terminalType Set the terminal type for the graphical output
     */
    DUNE_DEPRECATED_MSG("The signature of plot(string, int, string) has been changed to plot(string, string, int).")
    void plot(const std::string &title, const unsigned int plottingWindowNumber, const std::string& terminalType = "x11")
    {
        setTerminalType(terminalType);
        plot(title, terminalType, plottingWindowNumber);
    }

    /*!
     * \brief Plots the files for a specific window number, writes a gnuplot and png file.
     *
     * \param filename The name of the output file
     */
    void plot(const std::string &filename = "")
    {
        // set correct terminal and general options
        std::string plot = "reset\n";
        plot += "set datafile separator \'" + std::string(1, datafileSeparator_) + "\'\n";

        // set the terminal if the defaults were overwritten
        if (terminalType_.compare("x11") != 0 || linetype_.compare("solid") != 0)
          plot += "set term " + terminalType_ + " " + linetype_ + " " + " \n";

        // set the labels and axes ranges
        plot += "set xlabel \"" + xLabel_ + "\"\n";
        plot += "set ylabel \"" + yLabel_ + "\"\n";
        if (xRangeIsSet_)
            plot += "set xrange [" + std::to_string(xRangeMin_) + ":" + std::to_string(xRangeMax_) + "]" + "\n";
        if (yRangeIsSet_)
            plot += "set yrange [" + std::to_string(yRangeMin_) + ":" + std::to_string(yRangeMax_) + "]" + "\n";

        // set user defined options
        plot += plotOptions_ + "\n";

        // plot curves
        plot += "plot";
        for (unsigned int i = 0; i < curveFile_.size(); ++i)
        {
            plot += + " " + curveFile_[i]
                    + " " + curveOptions_[i]
                    + " title '" + curveTitle_[i] + "'";
            if (i < curveFile_.size()-1)
                plot += ",\\";
            plot += "\n";
        }

        // live plot of the results if gnuplot is installed
#ifdef HAVE_GNUPLOT
        if (openPlotWindow_)
            executeGnuplot(plot.c_str());
#endif

        // create a gnuplot file if a filename is specified
        if (filename.compare("") != 0)
        {
            plot += "\n";
            plot += "set term pngcairo size 800,600 " + linetype_ + " \n";
            plot += "set output \"" + filename + ".png\"\n";
            plot += "replot\n";
            std::string fileName = filename + ".gp";
            std::ofstream file;
            file.open(fileName);
            file << plot;
            file.close();
        }
    }

    /*!
     * \brief Restarts gnuplot
     */
    void resetAll(const bool persist = true)
    {
        close();
        open(persist);
        resetPlot();
    }

    /*!
     * \brief Deletes all plots from a plotting window and resets user-defined options
     */
    void resetPlot()
    {
        curveFile_.resize(0);
        curveOptions_.resize(0);
        curveTitle_.resize(0);
        plotOptions_ = "";
    }

    /*!
     * \brief Closes gnuplot
     */
    void open(const bool persist = true)
    {
        if (persist)
            pipe_ = popen((gnuplotPath_ + " -persist").c_str(), "w"); // "w" - writing
        else
            pipe_ = popen((gnuplotPath_).c_str(), "w");
    }

    /*!
     * \brief Closes gnuplot
     */
    void close()
    {
        if (pclose(pipe_) == -1)
            assert("Could not close pipe to Gnuplot!");
    }

    /*!
     * \brief Adds a function to list of plotted lines
     *
     * \param function Function to be plotted
     * \param plotName The name of the data set
     * \param plotOptions Specific gnuplot options passed to this plot
     */
    void addFunctionToPlot(const std::string function,
                           const std::string plotName,
                           const std::string plotOptions = "with lines")
    {
        curveFile_.push_back(function);
        curveOptions_.push_back(plotOptions);
        curveTitle_.push_back(plotName);
    }

    /*!
     * \brief Adds a file to list of plotted lines
     *
     * \param file Function to be plotted
     * \param plotName The name of the data set
     * \param plotOptions Specific gnuplot options passed to this plot
     */
    void addFileToPlot(const std::string file,
                       const std::string plotName,
                       const std::string plotOptions = "with lines")
    {
        curveFile_.push_back("'" + file + "'");
        curveOptions_.push_back(plotOptions);
        curveTitle_.push_back(plotName);
    }

    /*!
     * \brief Adds a data set and writes a data file
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
            checkNumber(x[i], "x[i] i=" + std::to_string(i) + " in " + fileName);
            checkNumber(y[i], "y[i] i=" + std::to_string(i) + " in " + fileName);
            file << x[i] << datafileSeparator_ << y[i] << std::endl;
        }
        file.close();

        // adding file to list of plotted lines
        curveFile_.push_back("'" + fileName + "'");
        curveOptions_.push_back(plotOptions);
        curveTitle_.push_back(plotName);
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
     * \param min The lowest plotted value for the x-axis
     * \param max The highest plotted value for the x-axis
     */
    void setXRange(Scalar min, Scalar max)
    {
        xRangeMin_ = min;
        xRangeMax_ = max;
        xRangeIsSet_ = true;
    }

    /*!
     * \brief Sets the range for the y-axis
     *
     * \param min The lowest plotted value for the y-axis
     * \param max The highest plotted value for the y-axis
     */
    void setYRange(Scalar min, Scalar max)
    {
        yRangeMin_ = min;
        yRangeMax_ = max;
        yRangeIsSet_ = true;
    }

    /*!
     * \brief Sets additional user-defined options
     *
     * \param option Additional line of option in gnuplot language
     */
    void setOption(std::string option)
    {
        plotOptions_ += option + "\n";
    }

    /*!
     * \brief Define whether the gnuplot window should be opened
     *
     * \param openPlotWindow Open gnuplot or not
     */
    void setOpenPlotWindow(bool openPlotWindow)
    {
        openPlotWindow_ = openPlotWindow;
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

    /*!
     * \brief Sets the terminal used for interactive outpu
     *
     * \param terminal The user-specified terminal
     */
    void setTerminalType(std::string terminal)
    {
        terminalType_ = terminal;
    }

    /*!
     * \brief Use dashed (true) or solid (false) lines
     *
     * \param dashed Use dashed lines
     */
    void useDashedLines(bool dashed)
    {
        linetype_ = dashed ? "dashed" : "solid";
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

    std::FILE * pipe_;
    bool openPlotWindow_;
    bool persist_;
    std::string terminalType_;
    char datafileSeparator_;
    std::string linetype_;
    StringVector curveFile_;
    StringVector curveOptions_;
    StringVector curveTitle_;
    bool interaction_;
    Scalar xRangeMin_;
    Scalar xRangeMax_;
    Scalar xRangeIsSet_;
    Scalar yRangeMin_;
    Scalar yRangeMax_;
    Scalar yRangeIsSet_;
    std::string xLabel_;
    std::string yLabel_;
    std::string plotOptions_;
    std::string gnuplotPath_;
};
} // end of namespace
#endif // DUMUX_GNUPLOT_INTERFACE_HH
