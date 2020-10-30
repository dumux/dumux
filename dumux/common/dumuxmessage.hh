// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Common
 * \brief Provides the class creating the famous DuMu<sup>x</sup> start and end messages
 */
#ifndef DUMUX_MESSAGE_HH
#define DUMUX_MESSAGE_HH

#include <iomanip>
#include <iostream>
#include <ctime>

#include <dumux/io/format.hh>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief DuMu<sup>x</sup> start and end message.
 */
class DumuxMessage
{
    //! The current number of messages. Please adjust if you add one.
    static const int nMessages_ = 12;

public:

    /*!
     * \brief Selects random messages to write out at the start and end of a simulation run.
     * \param firstCall Indicates if it's the first call and we have to dice (simulation is starting).
     */
    static void print(bool firstCall = false)
    {
        // initialize in case someone forgets to set first call
        static int dice = 8;

        if(firstCall)
        {
            // roll the dice to decide which start message will be displayed
            std::srand(std::time(0));
            dice = std::rand() % (nMessages_ + 1);
        }

        std::cout << std::endl;

        switch (dice)
        {
            case 0:
                if(firstCall)
                    std::cout << "Welcome aboard DuMuX airlines. Please fasten your seatbelts! "
                              << "Emergency exits are near the time integration." << std::endl;
                else
                    std::cout << "We hope that you enjoyed simulating with us " << std::endl
                              << "and that you will choose us next time, too." << std::endl;
            break;
            case 1:
                if(firstCall)
                    std::cout << "Let's get the cow off the ice." << std::endl;
                else
                    std::cout << "DuMuX got the cow off the ice." << std::endl;
            break;
            case 2:
                if(firstCall)
                    std::cout << "Science, my lad, is made up of mistakes, but they are "
                              << "mistakes which it is useful to make, because they lead little "
                              << "by little to the truth." << std::endl
                              << " - Jules Verne, A journey to the center of the earth" << std::endl;
                else
                    std::cout << "[We see that] science is eminently perfectible, and that each theory has "
                              << "constantly to give way to a fresh one." << std::endl
                              << " - Jules Verne, Journey to the Center of the Earth" << std::endl;

            break;
            case 3:
                if(firstCall)
                    std::cout << "Wherever he saw a hole he always wanted to know the depth of it. "
                              << "To him this was important." << std::endl
                              << " - Jules Verne, A journey to the center of the earth" << std::endl;
                else
                    std::cout << "We may brave human laws, but we cannot resist natural ones." << std::endl
                              << " - Jules Verne, 20,000 Leagues Under the Sea" << std::endl;
            break;
            case 4:
                if(firstCall)
                    std::cout << "Silence - to delight Bernd." << std::endl;
                else
                    std::cout << std::endl << std::endl;
            break;
            case 5:
                    std::cout << "Don't panic... !" << std::endl;
            break;
            case 6:
                if(firstCall)
                    std::cout << "You idiot! You signed the order to destroy Earth!" << std::endl
                              << " - Douglas Adams, HGttG" << std::endl;
                else
                    std::cout << "Marvin: I've been talking to the main computer." << std::endl
                              << "Arthur: And?" << std::endl
                              << "Marvin: It hates me." << std::endl
                              << " - Douglas Adams, HGttG" << std::endl;
            break;
            case 7:
                if(firstCall)
                    std::cout << "In the beginning the Universe was created. This has made a lot of "
                              << "people very angry and has been widely regarded as a bad move.!" << std::endl
                              << " - Douglas Adams, HGttG " << std::endl;
                else
                    std::cout << "Forty-two. I checked it very thoroughly, and that quite definitely is the answer. I think "
                              << "the problem, to be quite honest with you, is that you\'ve never actually known what the question is." << std::endl
                              << " - Douglas Adams, HGttG " << std::endl;
            break;
            case 8:
                std::cout << "                ##                  @@@@          @   @         @ @" << std::endl;
                std::cout << "             ###   #                @   @         @@ @@          @ " << std::endl;
                std::cout << "           ##       #               @   @  @   @  @ @ @  @   @  @ @" << std::endl;
                std::cout << "         ##          #              @   @  @   @  @   @  @   @     " << std::endl;
                std::cout << "        #             #             @@@@    @@@   @   @   @@@      " << std::endl;
                std::cout << "       #               #                                           " << std::endl;
                std::cout << "      #                 #                                          " << std::endl;
                std::cout << "     #                   ##        %%%                     " << std::setw(8) << std::right << DUMUX_VERSION << std::endl;
                std::cout << "    #                      ###    %   %  %%     %%                 " << std::endl;
                std::cout << "####                          #%%%     %%  %%%%%  %%%%%%%%%%%%%%%%%" << std::endl;
            break;
            case 9:
                std::cout << "###         #   #        # #                            " << std::endl;
                std::cout << "#  #  #  #  ## ##  #  #   #                             " << std::endl;
                std::cout << "#  #  #  #  # # #  #  #  # #                            " << std::endl;
                std::cout << "###    ##   #   #   ##                                  " << std::endl;
                std::cout << "                                                        " << std::endl;
                std::cout << "Dune for Multi-{ Phase,                                 " << std::endl;
                std::cout << "                 Component,                             " << std::endl;
                std::cout << "                 Scale,                                 " << std::endl;
                std::cout << "                 Physics,                               " << std::endl;
                std::cout << "                 ...} flow and transport in porous media" << std::endl;
            break;
            case 10:
                if(firstCall)
                    std::cout << "Elliot Carver: Mr. Jones, are we ready to release our new software?" << std::endl
                              << "Jones: Yes, sir. As requested, it's full of bugs, which means people will be forced to upgrade for years." << std::endl
                              << " - James Bond, Tomorrow Never Dies" << std::endl;
                else
                {
                    std::cout << "Elliot Carver: Outstanding." << std::endl
                              << " - James Bond, Tomorrow Never Dies" << std::endl;
                }
            break;
            case 11:
                if(firstCall)
                    std::cout << "Chuck Norris has successfully compiled DuMuX." << std::endl;
                else
                    std::cout << "Chuck Norris has compiled DuMuX even two times in a row!" << std::endl;
            break;
            case 12:
                if (firstCall)
                    std::cout << Fmt::format("  ┌{0:─^{2}}┐\n"
                                             "  │{1: ^{2}}│\n"
                                             "  └{0:─^{2}}┘\n", "", Fmt::format("DuMuX {} \u2661", DUMUX_VERSION), 20);
                else
                    std::cout << "\n" << std::endl;
            break;

            // Note: If you add a case, you have to increase the number of messages (nMessages_ variable).

            default:    // silence to delight Bernd
                return;
        }
        std::cout << std::endl;
    }
};

} // end namespace Dumux

#endif
