import os
import sys
import functools
import subprocess
import traceback


def getCommandErrorHints(command):
    if "git " in command:
        return "It seems that a git command failed. Please check:\n" \
               "    -- is the module registered as git repository?\n" \
               "    -- is upstream defined for the branch?"
    return None


# execute a command and retrieve the output
def runCommand(command, suppressTraceBack=False, errorMessage=''):
    try:
        return subprocess.run(command,
                              shell=True, check=True,
                              text=True, capture_output=True).stdout
    except Exception:
        eType, eValue, eTraceback = sys.exc_info()
        if suppressTraceBack:
            traceback.print_exception(eType, eType(errorMessage), None)
        elif errorMessage:
            traceback.print_exception(eType, eType(errorMessage), eTraceback)
        else:
            print("An error occurred during subprocess run:")
            print("-- command: {}".format(command))
            print("-- folder: {}".format(os.getcwd()))
            traceback.print_exception(eType, eValue, eTraceback)
            hints = getCommandErrorHints(command)
            if hints is not None:
                print(hints)


# decorator to call function from within the given path
def callFromPath(path):
    def decorator_callFromPath(callFunc):
        @functools.wraps(callFunc)
        def wrapper_callFromPath(*args, **kwargs):
            curPath = os.getcwd()
            os.chdir(path)
            result = callFunc(*args, **kwargs)
            os.chdir(curPath)
            return result
        return wrapper_callFromPath
    return decorator_callFromPath
