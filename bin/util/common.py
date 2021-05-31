import os
import sys
import functools
import subprocess


# execute a command and retrieve the output
def runCommand(command, suppressTraceBack=False, errorMessage=''):
    try:
        return subprocess.run(command,
                              shell=True, check=True,
                              text=True, capture_output=True).stdout
    except Exception:
        if suppressTraceBack:
            sys.excepthook(Exception, Exception(errorMessage), None)
        elif not errorMessage:
            print("An error occurred during subprocess run:")
            print("-- command: {}".format(command))
            print("-- folder: {}".format(os.getcwd()))
            print("-- error: {}".format(sys.exc_info()[1]))
            if "git " in command:
                print()
                print("It seems that a git command failed. Please check:\n"
                      "    -- is the module registered as git repository?\n"
                      "    -- is upstream defined for the branch?\n")
            raise
        else:
            raise Exception(errorMessage)


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
