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
 * \ingroup InputOutput
 * \author Timo Koch
 * \brief Tools for creating logs
 */
#ifndef DUMUX_IO_LOG_HH
#define DUMUX_IO_LOG_HH

#include <iostream>
#include <string>
#include <type_traits>

#include <dune/common/parallel/mpihelper.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

//! log levels in Dumux
enum class LogLevel : int
{ warning = 40, info = 30, progress = 20, trace = 10, debug = 0 };

class Logger
{
public:
    // create a new logger
    Logger()
    : indentationLevel_(0)
    , ostream_(&std::cout)
    , minLogLevel_(LogLevel::progress)
    , sleeping_(false)
    , comm_(Dune::MPIHelper::getCollectiveCommunication())
    {}

    /*
     * \brief Log a new message
     */
    void log(const std::string& message, LogLevel logLevel)
    {
        assert (static_cast<int>(logLevel) < static_cast<int>(LogLevel::warning));
        writeLog_(message, logLevel);
    }

    /*
     * \brief Log a warning
     */
    void warning(const std::string& warning)
    {
        const std::string message = ">>> Dumux warning: " + warning;
        writeLog_(message, LogLevel::warning);
    }

    /*
     * \brief Begin a new task and inrease indentation level
     * \note Do not forget to call endTask() as well
     */
    void beginTask(const std::string& taskDescription, LogLevel logLevel)
    {
        if (logLevelApplies_(logLevel))
        {
            writeLog_(taskDescription, logLevel);
            ++indentationLevel_;
        }
    }

    /*
     * \brief End the current task and decrease indentation level
     * \note Only call this after calling beginTask()
     */
    void endTask()
    {
        --indentationLevel_;
    }

    /*
     * \brief Disable the logger temporarily
     */
    void setSleeping(bool sleeping = true)
    { sleeping_ = sleeping; }

    /*
     * \brief Return if the logger is currently sleeping (inactive)
     */
    bool isSleeping() const
    { return sleeping_; }

    /*
     * \brief Set the current output stream
     */
    void setOutputStream(std::ostream& ostream)
    { ostream_ = &ostream; }

    /*
     * \brief Set the minimum log level
     */
    void setLogLevel(LogLevel logLevel)
    { minLogLevel_ = logLevel; }

private:
    /*
     * \brief if a given log level has a high enough priority
     */
    bool logLevelApplies_(LogLevel logLevel) const
    { return (static_cast<int>(logLevel) >= static_cast<int>(minLogLevel_)); }

    /*
     * \brief write a log entry
     */
    void writeLog_(std::string message, LogLevel logLevel) const
    {
        // if we are sleeping or the log level is too low skip
        if (sleeping_ || !logLevelApplies_(logLevel))
            return;

        // only write out errors and worse on all processes
        if (comm_.rank() > 0 && static_cast<int>(logLevel) < static_cast<int>(LogLevel::warning))
            return;

        // in parallel prepend the processor rank
        if (comm_.size() > 1 && static_cast<int>(logLevel) >= static_cast<int>(LogLevel::warning))
            message = "Process " + std::to_string(comm_.rank()) + ": " + message;

        // set the indentation level
        std::stringstream indent;
        for (int i = 0; i < indentationLevel_; ++i)
            indent << "  ";
        message = indent.str() + message;

        // write out and new line and flush
        *ostream_ << message << std::endl;
    }

    int indentationLevel_;
    std::ostream* ostream_;
    LogLevel minLogLevel_;
    bool sleeping_;

    Dune::CollectiveCommunication<typename Dune::MPIHelper::MPICommunicator> comm_;
};

//! The log manager that manages a logger singleton
struct LogManager
{
    LogManager() = delete;
    LogManager(const LogManager&) = delete;
    void operator=(LogManager const&) = delete;

    static Logger& logger()
    {
        static Logger logger;
        return logger;
    }
};

//! The log namespace containing all convenience functions
namespace Log {

struct LogTag {};
struct Debug : public LogTag { static constexpr LogLevel logLevel = LogLevel::debug; };
struct Trace : public LogTag { static constexpr LogLevel logLevel = LogLevel::trace; };
struct Progress : public LogTag { static constexpr LogLevel logLevel = LogLevel::progress; };
struct Info : public LogTag { static constexpr LogLevel logLevel = LogLevel::info; };
struct Warning : public LogTag { static constexpr LogLevel logLevel = LogLevel::warning; };
struct End {};

//! helper to collect the whole log message before passing it to the logger
struct LogCollector
{
    explicit LogCollector(LogLevel l) noexcept
    : logLevel(l)
    {}

    std::stringstream stream;
    LogLevel logLevel;
};

//! helper to collect the whole log message before passing it to the logger
template<class Input>
LogCollector&& operator<< (LogCollector&& collector, const Input& input)
{
    collector.stream << input;
    return std::move(collector);
}

//! helper to collect the whole log message before passing it to the logger
void operator<< (LogCollector&& collector, const End& end)
{
    if (LogManager::logger().isSleeping())
        return;

    if (collector.logLevel == LogLevel::warning)
        LogManager::logger().warning(collector.stream.str());
    else
        LogManager::logger().log(collector.stream.str(), collector.logLevel);
}

constexpr Debug debug;
constexpr Trace trace;
constexpr Progress progress;
constexpr Info info;
constexpr Warning warning;
constexpr End end;

/*
 * \brief Attach logs to the logger
 * Examples:
 * (1)  Log::warning << "Computing the composition without a constraint solver object is deprecated!" << Log::end;
 * (2)  Log::trace << "Found " << 3 << " intersections!" << Log::end;
 */
template<class Tag, class Input, typename std::enable_if_t<std::is_base_of<LogTag, Tag>::value, int> = 0>
LogCollector operator<< (const Tag& tag, const Input& input)
{
    LogCollector collector(tag.logLevel);
    collector.stream << input;
    return collector;
}

/*
 * \brief Begin a new task and inrease indentation level
 * \note Do not forget to call endTask() as well
 */
inline void beginTask(const std::string& taskDescription, LogLevel logLevel)
{
    if (LogManager::logger().isSleeping())
        return;

    LogManager::logger().beginTask(taskDescription, logLevel);
}

/*
 * \brief End the current task and decrease indentation level
 * \note Only call this after calling beginTask()
 */
inline void endTask()
{
    if (LogManager::logger().isSleeping())
        return;

    LogManager::logger().endTask();
}

/*
 * \brief Disable the logger temporarily
 */
inline void setSleeping(bool sleeping = true)
{ LogManager::logger().setSleeping(sleeping); }

/*
 * \brief Return if the logger is currently sleeping (inactive)
 */
inline bool isSleeping()
{ return LogManager::logger().isSleeping(); }

/*
 * \brief Set the current output stream
 */
inline void setOutputStream(std::ostream& ostream)
{ LogManager::logger().setOutputStream(ostream); }

} // end namespace Log

/*
 * \brief A message that gets printed multiple times
 * Examples:
 * LogEventMessage intersectionMsg("Found intersection");
 * for (const auto& e : elements(gridView))
 *     if (foundIntersection(e))
 *         intersectionMsg();
 */
class LogEventMessage
{
public:
    LogEventMessage(const std::string& baseMsg, LogLevel logLevel, std::size_t maxCount = 1)
    : count_(0)
    , maxCount_(maxCount)
    , baseMsg_(baseMsg)
    , logLevel_(logLevel)
    {}

    void operator() (const std::string& add = "")
    {
        if (count_ < maxCount_)
            LogManager::logger().log(baseMsg_ + " " + add, logLevel_);

        ++count_;
    }

    std::size_t count() const
    { return count_; }

    void summarize() const
    {
        LogManager::logger().log("Event \"" + baseMsg_ + "\" occured " + std::to_string(count()) + " times.", logLevel_);
    }

private:
    std::size_t count_;
    std::size_t maxCount_;
    std::string baseMsg_;
    LogLevel logLevel_;
};

} // end namespace Dumux

#endif
