// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Evaluating string math expressions containing named variables
 */
#ifndef DUMUX_COMMON_FUNCTION_FROM_STRING_EXPRESSION_HH
#define DUMUX_COMMON_FUNCTION_FROM_STRING_EXPRESSION_HH

#include <array>
#include <mutex>
#include <string>
#include <string_view>
#include <type_traits>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/io/format.hh>
#include <dumux/io/expression/exprtk.hpp>

namespace Dumux {

/*!
 * \ingroup Core
 * \brief Evaluating string math expressions containing named variables
 * \tparam numVars number of variables in the expression; number of function arguments of the call operator
 * \tparam Scalar type of numerical values in the expression
 *
 * Example usage
 \code{.cpp}
    // Create a callable f(x,t) from a string expression containing variable literals.
    // The constructor compiles the expression in the constructor making calls efficient.
    // The constructor throws a Dune::IOError with detailed info if parsing fails.
    std::string expr = getParam("Problem.Function"); // e.g. "5*x + x*sin(x*t)"
    FunctionFromStringExpression<2> f(expr, "xt"); // variables "x" and "t"

    // evaluate function, result is double (the default scalar type)
    const double x = 1.0, t = 2.0;
    const double result = f(x, t);
 \endcode
 *
 * For variables with several characters construct
 \code{.cpp}
    // variables "pos" and "time"
    FunctionFromStringExpression<2> f(expr, std::array<std::string, 2>{{"pos", "time"}});
 \endcode
 *
 */
template<std::size_t numVars, class Scalar = double>
class FunctionFromStringExpression
{
    using SymbolTable = exprtk::symbol_table<Scalar>;
    using Expression = exprtk::expression<Scalar>;
    using Parser = exprtk::parser<Scalar>;

public:
    static constexpr std::size_t numVariables = numVars;

    //! \brief Constructor from math expression and array of variable names
    FunctionFromStringExpression(const std::string& expression, const std::array<std::string, numVars>& variableNames)
    { initialize_(expression, variableNames); }

    //! \brief Delegating constructor using all characters of a string as variables
    //! \note Calling FunctionFromStringExpression(expr, "xt") uses "x" and "t" as variable names
    FunctionFromStringExpression(const std::string& expression, std::string_view variableNames)
    : FunctionFromStringExpression(expression, extractVariableNames_(variableNames, std::make_index_sequence<numVars>{})) {}

    template<class S, std::enable_if_t<std::is_convertible_v<Scalar, S>, int> = 0>
    Scalar operator() (const std::array<S, numVars>& params) const
    { return evalRandomAcessImpl_(params); }

    template<class S, std::enable_if_t<std::is_convertible_v<Scalar, S>, int> = 0>
    Scalar operator() (const Dune::FieldVector<S, numVars>& params) const
    { return evalRandomAcessImpl_(params); }

    template<class ...Params, std::enable_if_t<(sizeof...(Params) == numVars) && (std::is_convertible_v<Scalar, std::decay_t<Params>> && ...), int> = 0>
    Scalar operator() (Params&&... params) const
    { return evalRandomAcessImpl_(std::array<Scalar, numVars>{ std::forward<Params>(params)... }); }

    void setVerbosity(unsigned int v)
    { verbosity_ = v; }

private:
    template<class RandomAccessContainer>
    Scalar evalRandomAcessImpl_(const RandomAccessContainer& params) const
    {
        std::lock_guard lock(evalMutex_);
        for (std::size_t i = 0; i < numVars; ++i)
            variables_[i] = params[i];
        return expression_.value();
    }

    template<std::size_t... I>
    std::array<std::string, numVars> extractVariableNames_(std::string_view names, std::index_sequence<I...>) const
    {
        static_assert(numVars == sizeof...(I), "Number of variables has to match size of index set.");
        if (names.size() != numVars)
            DUNE_THROW(Dune::IOError, "Number of variables in '"
                << names << "' does not match number of function arguments: " << numVars);
        return { std::string(1, names.at(I))... };
    }

    //! Parse the math expression and throw detailed error message when this fails
    void initialize_(const std::string& expression, const std::array<std::string, numVars>& variableNames)
    {
        for (std::size_t i = 0; i < numVars; ++i)
            symbolTable_.add_variable(std::string{variableNames[i]}, variables_[i]);
        symbolTable_.add_constants();
        expression_.register_symbol_table(symbolTable_);

        if (!parser_.compile(expression, expression_))
        {
            std::stringstream ss;
            ss << Fmt::format("Parsing expression '{}' failed.\n", expression);

            if (verbosity_ >= 1)
            {
                for (std::size_t i = 0; i < parser_.error_count(); ++i)
                {
                    const auto error = parser_.get_error(i);

                    ss << Fmt::format(
                        "-- error (position: {:02d}, type: {}): {}\n",
                        static_cast<unsigned int>(error.token.position),
                        exprtk::parser_error::to_str(error.mode).c_str(),
                        error.diagnostic.c_str()
                    );
                }
            }

            DUNE_THROW(Dune::IOError, ss.str());
        }
        else if (verbosity_ >= 2)
        {
            std::cout << Fmt::format(
                "Successfully parsed math expression '{}'\n",
                expression
            );
        }
    }

    unsigned int verbosity_ = 2;
    SymbolTable symbolTable_;
    Expression expression_;
    Parser parser_;
    mutable std::array<Scalar, numVars> variables_;
    mutable std::mutex evalMutex_;
};

} // end namespace Dumux

#endif
