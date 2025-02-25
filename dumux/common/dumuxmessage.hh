// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
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
 * \ingroup Core
 * \brief DuMu<sup>x</sup> start and end message.
 */
class DumuxMessage
{
    //! The current number of messages. Please adjust if you add one.
    static const int nMessages_ = 34;

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
                {
                    std::cout << "  ┌──────────────────┐" << std::endl;
                    std::cout << Fmt::format("  │{:^20}│", Fmt::format("DuMuX {} \u2661", DUMUX_VERSION)) << std::endl;
                    std::cout << "  └──────────────────┘" << std::endl;
                }
                else
                    std::cout << "\n" << std::endl;
            break;
            case 13:
                if(firstCall)
                {
                    std::cout << "Everything starts somewhere, though many physicists disagree." << std::endl
                              << " - Terry Pratchett " << std::endl;
                }
                else
                {
                    std::cout << "Opera happens because a large number of things amazingly fail to go wrong." << std::endl
                              << " - Terry Pratchett " << std::endl;
                }
            break;
            case 14:
                    std::cout << "To infinity and beyond." << std::endl
                              << " - Buzz Lightyear, Toy Story" << std::endl;
            break;
            case 15:
                if(firstCall)
                {
                    std::cout << "C makes it easy to shoot yourself in the foot; C++ makes it harder, but when you do it blows your whole leg off." << std::endl
                              << " - Bjarne Stroustrup " << std::endl;
                }
                else
                {
                    std::cout << "There's an old story about the person who wished his computer were as easy to use as his telephone." << std::endl
                              << "That wish has come true, since I no longer know how to use my telephone." << std::endl
                              << " - Bjarne Stroustrup " << std::endl;
                }
            break;
            case 16:
                if(firstCall)
                {
                    std::cout << "Now, all we need is a little Energon and a lot of luck!" << std::endl
                              << " - Optimus Prime, The Transformers: The Movie " << std::endl;
                }
                else
                {
                    std::cout << "Sometimes even the wisest of men and machines can be in error." << std::endl
                              << " - Optimus Prime, The Transformers: The Movie " << std::endl;
                }
            break;
            case 17:
                if(firstCall)
                {
                    std::cout << "Let's go. In and out, 20 minutes adventure." << std::endl
                              << " - Rick Sanchez, Rick & Morty " << std::endl;
                }
                else
                {
                    std::cout << "Losers look stuff up while the rest of us are carpin' all them diems." << std::endl
                              << " - Summer Smith, Rick & Morty" << std::endl;
                }
            break;
            case 18:
                if(firstCall)
                {
                    std::cout << "It's the job that's never started as takes longest to finish." << std::endl
                              << " - Sam Gamgee, LotR " << std::endl;
                }
                else
                {
                    std::cout << "He that breaks a thing to find out what it is, has left the path of wisdom." << std::endl
                              << " - Gandalf, LotR " << std::endl;
                }
            break;
            case 19:
                if(firstCall)
                {
                    std::cout << "The Ring has awoken, it's heard its master's call." << std::endl
                              << " - Gandalf, LotR " << std::endl;
                }
                else
                {
                    std::cout << "It's a dangerous business, Frodo, going out your door. " << std::endl
                              << "You step onto the road, and if you don't keep your feet, there's no knowing where you might be swept off to." << std::endl
                              << " - Frodo Baggins, LotR " << std::endl;
                }
            break;
            case 20:
                if(firstCall)
                {
                    std::cout << "Who knows? Have patience. Go where you must go, and hope!" << std::endl
                              << " - Gandalf, LotR " << std::endl;
                }
                else
                {
                    std::cout << "Don't adventures ever have an end? I suppose not. Someone else always has to carry on the story." << std::endl
                              << " - Bilbo Baggins, LotR " << std::endl;
                }
            break;
            case 21:
                if(firstCall)
                {
                    std::cout << "As long as I'm better than everyone else I suppose it doesn't matter." << std::endl
                              << " - Jamie Lannister, GoT" << std::endl;
                }
                else
                {
                    std::cout << "My watch has ended." << std::endl
                              << " - Jon Snow, GoT" << std::endl;
                }
            break;
            case 22:
                if(firstCall)
                {
                    std::cout << "You'll find I'm full of surprises." << std::endl
                              << " - Luke Skywalker, Star Wars: The Empire Strikes Back " << std::endl;
                }
                else
                {
                    std::cout << "I find your lack of faith disturbing." << std::endl
                              << " - Darth Vader, Star Wars: A New Hope " << std::endl;
                }
            break;
            case 23:
                if(firstCall)
                {
                    std::cout << "Here goes nothing." << std::endl
                              << " - Lando Calrissian, Star Wars: Return of the Jedi" << std::endl;
                }
                else
                {
                    std::cout << "Chewie, we're home." << std::endl
                              << " - Han Solo, Star Wars: The Force Awakens" << std::endl;
                }
            break;
            case 24:
                if(firstCall)
                {
                    std::cout << "The Force is strong with this one." << std::endl
                              << " - Darth Vader, Star Wars: A New Hope " << std::endl;
                }
                else
                {
                    std::cout << "In my experience, there's no such thing as luck." << std::endl
                              << " - Obi-Wan Kenobi, Star Wars: A New Hope " << std::endl;
                }
            break;
            case 25:
                if(firstCall)
                {
                    std::cout << "The city's central computer told you? R2D2, you know better than to trust a strange computer!" << std::endl
                              << " - C3PO, Star Wars: The Empire Strikes Back " << std::endl;
                }
                else
                {
                    std::cout << "He's quite clever, you know...for a human being." << std::endl
                              << " - C3PO, Star Wars: The Empire Strikes Back " << std::endl;
                }
            break;
            case 26:
                if(firstCall)
                {
                    std::cout << "I know some things. I can, you know, do math and stuff." << std::endl
                              << " - Harry Potter " << std::endl;
                }
                else
                {
                    std::cout << "Harry then did something that was both very brave and very stupid." << std::endl
                              << " - Harry Potter and the Sorcerer's Stone " << std::endl;
                }
            break;
            case 27:
                if(firstCall)
                {
                    std::cout << "I'll be in my bedroom, making no noise and pretending I'm not there." << std::endl
                              << " - Harry Potter " << std::endl;
                }
                else
                {
                    std::cout << "Honestly, if you were any slower, you'd be going backward." << std::endl
                              << " - Draco Malfoy " << std::endl;
                }
            break;
            case 28:
                    std::cout << "I can do this all day." << std::endl
                              << " - Captain America " << std::endl;
            break;
            case 29:
                if(firstCall)
                {
                    std::cout << "Your scientists were so preoccupied with whether or not they could, they didn't stop to think if they should." << std::endl
                              << " - Ian Malcolm, Jurassic Park " << std::endl;
                }
                else
                {
                    std::cout << "Boy, do I hate being right all the time." << std::endl
                              << " - Ian Malcolm, Jurassic Park " << std::endl;
                }
            break;
            case 30:
                if(firstCall)
                {
                    std::cout << "It's a UNIX System! I know this! "
                              << " - Lex Murphy, Jurassic Park " << std::endl;
                }
                else
                {
                    std::cout << "When you gotta go, you gotta go." << std::endl
                              << " - Ian Malcolm, Jurassic Park " << std::endl;
                }
            break;
            case 31:
                if(firstCall)
                {
                    std::cout << "Whatever happens, that's the plan. "
                              << " - Kayla Watts, Jurassic World Dominion " << std::endl;
                }
                else
                {
                    std::cout << "Can we start over?" << std::endl
                              << " - Claire Dearing, Jurassic World Dominion " << std::endl;
                }
            break;
            case 32:
                if(firstCall)
                {
                    std::cout << "The code is more what you'd call 'guidelines' than actual rules. "
                              << " - Hector Barbossa, Pirates of the Caribbean " << std::endl;
                }
                else
                {
                    std::cout << "Did everyone see that? Because I will not be doing it again." << std::endl
                              << " - Jack Sparrow, Pirates of the Caribbean " << std::endl;
                }
            break;
            case 33:
                if(firstCall)
                {
                    std::cout << "If you were waiting for the opportune moment, that was it. "
                              << " - Jack Sparrow, Pirates of the Caribbean " << std::endl;
                }
                else
                {
                    std::cout << "I love those moments. I like to wave at them as they pass by." << std::endl
                              << " - Jack Sparrow, Pirates of the Caribbean " << std::endl;
                }
            break;
            case 34:
                std::cout << "And that was without even a single drop of rum." << std::endl
                          << " - Jack Sparrow, Pirates of the Caribbean " << std::endl;
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
