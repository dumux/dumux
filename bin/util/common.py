import os
import re
import sys
import fnmatch
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


# query a yes/no answer from the user
def query_yes_no(question, default="yes"):
    affirmative = ["yes", "y", "ye"]
    negative = ["no", "n"]

    def get_choices():
        return ", ".join(c for c in affirmative + negative)

    def is_affirmative(choice): return choice in affirmative
    def is_negative(choice): return choice in negative
    def is_valid(choice): return is_affirmative(choice) or is_negative(choice)

    if not is_valid(default):
        raise ValueError("\nInvalid default answer: '{}', choices: '{}'\n"
                         .format(default, get_choices()))

    if default is None:
        prompt = " [y/n] "
    else:
        prompt = " [Y/n] " if is_affirmative(default) else " [y/N] "

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()

        if default is not None and choice == "":
            return True if is_affirmative(default) else False

        if not is_valid(choice):
            sys.stdout.write("\nInvalid answer: '{}'. Choose from '{}'.\n"
                             .format(choice, get_choices()))

        return True if is_affirmative(choice) else False


def header_filter_cpp():
    return lambda file_name: file_name == 'config.h'


# get all project headers included by a cpp file
def get_included_project_headers_cpp(file,
                                     project_base,
                                     headers=[],
                                     header_filter=header_filter_cpp()):
    file_path = os.path.join(project_base, file)
    if not os.path.exists(file_path):
        raise IOError(f'Cpp file {file_path} does not exist')

    with open(file_path, 'r') as f:
        content = f.read()
        header_in_bracket = re.findall(r'#include\s+<(.+?)>', content)
        header_in_quotation = re.findall(r'#include\s+"(.+?)"', content)

        def process(path_in_project):
            header_path = os.path.join(project_base, path_in_project)
            if os.path.exists(header_path):
                if not header_filter(path_in_project):
                    if header_path not in headers:
                        headers.append(header_path)
                        get_included_project_headers_cpp(
                            header_path, project_base,
                            headers, header_filter
                        )

        for header in header_in_bracket:
            process(header)
        for header in header_in_quotation:
            abs_header_path = os.path.join(os.path.dirname(file), header)
            project_path = os.path.relpath(abs_header_path, project_base)
            process(project_path)
    return headers


# find all files below the given folder that match the given pattern
def findMatchingFiles(path, pattern):
    result = []
    for root, dirs, files in os.walk(path):
        root_rel = os.path.relpath(root, path)
        for file in files:
            if fnmatch.fnmatch(file, pattern):
                result.append(os.path.join(root_rel, file))
    return result
