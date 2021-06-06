#!/usr/bin/env python3
"""
This script extracts some specified applications into a separate Dune module.
For example make a dumux-pub repository accompanying a scientific paper.
"""
import sys
import os
import subprocess
import os.path
import argparse
import shutil
from distutils.dir_util import copy_tree
import re
import multiprocessing as mp
import fnmatch
import itertools
from functools import partial
from util import getPersistentVersions
from util import versionTable
from makeinstallscript import makeInstallScript
try:
    path = os.path.split(os.path.abspath(__file__))[0]
    sys.path.append(os.path.join(path, '../bin/util'))
    from getmoduleinfo import getDependencies
    from common import callFromPath, runCommand
except Exception:
    sys.exit('Could not import common modul or getModuleInfo')


# return the list of included headers including the header itself
def add_headers_recursively(header_path, curret_headers, module_path):
    hh = []
    if os.path.exists(header_path):
        if header_path not in curret_headers:
            hh.append(header_path)
            hh += search_headers(header_path, module_path)
    return hh


# function to search matching header files
def search_headers(source_file, module_path):
    headers = []

    with open(source_file, 'r') as f:
        content = f.read()
        header_in_bracket = re.findall(r'#include\s+<(.+?)>', content)
        header_in_quotation = re.findall(r'#include\s+"(.+?)"', content)

    # search for includes relative to the module path
    for header in header_in_bracket + header_in_quotation:
        header_path = os.path.join(module_path, header)
        headers += add_headers_recursively(header_path, headers, module_path)

    # search for includes relative to the path of the including file
    # only allow quoted includes for this
    for header in header_in_quotation:
        if header == "config.h":
            continue
        header_path = os.path.join(os.path.dirname(source_file), header)
        headers += add_headers_recursively(header_path, headers, module_path)

    return headers


# function asking user to answer yes or no question
def query_yes_no(question, default="yes"):
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)
    while True:
        sys.stdout.write("\n\n" + question + prompt)
        choice = input().lower()
        if default is not None and choice == "":
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


# query for git remote url and make sure it is not empty
def get_remote_url(repo_path):
    run_from_mod = callFromPath(repo_path)(runCommand)
    while True:
        if query_yes_no("Do you own a subfolder in dumux-pub?"):
            nameyearx = input("Type below the information of"
                              " AuthorLastNameYearx (e.g. Luigi2020b)\n")
            remoteurl = 'https://git.iws.uni-stuttgart.de'\
                        '/dumux-pub/{}.git'.format(nameyearx.lower())
        else:
            remoteurl = input("Provide URL of your remote repository:\n")
        check_remote_repo = run_from_mod('git ls-remote {}'.format(remoteurl))
        if (check_remote_repo is None):
            sys.stdout.write("\nERROR: Please re-enter correct information.\n")
        elif (check_remote_repo == ''):
            return remoteurl
        else:
            sys.stdout.write("ERROR: The remote reposity is not empty!.\n")


def extract_sources_files(module_dir, subfolders):
    def find_files_with_pattern(pattern, path):
        result = []
        for root, dirs, files in os.walk(path):
            for name in files:
                if fnmatch.fnmatch(name, pattern):
                    result.append(os.path.join(root, name))
        return result

    sources = []
    for folder in subfolders:
        folder_path = os.path.join(module_dir, folder)
        if not os.path.isdir(folder_path):
            raise NameError(f"Subfolder '{folder}' is not a subfolder of '{module_dir}'")
        sources += find_files_with_pattern("*.cc", os.path.abspath(folder_path))
    return sources


def check_module(module_name):
    try:
        with open(f"{module_name}/dune.module") as module_file:
            content = module_file.read()
            if f"Module: {module_name}" not in content:
                raise Exception(f"Invalid dune.module in {module_name}")
    except OSError:
        print("Could not find new Dune module. Aborting")
        raise


###################################################################
# Some general information for users of the script
###################################################################
def info_explanations(module_dir, module_path, subfolders, source_files):
    return f"""
This script will extract the following subfolders of
the module '{module_dir}':
{os.linesep.join([f"   - {os.path.relpath(f, module_path)}" for f in subfolders])}

and all headers contained in '{module_dir}'
tha are required to build the exectutables from the sources:
{os.linesep.join([f"   - {s}" for s in source_files])}

The extracted files are copied into a new DUNE module
retaining the directory structure. Required files for
creating a working DUNE module (like CMakeLists.txt)
will be created and/or updated.

The script 'duneproject' will run now. Please answer all
upcoming queries to the best of your knowledge.

Important: the new module should NOT depend on the module '{module_dir}'
"""


###################################################################
# Main part of README.md
###################################################################
def info_readme_main(module_dir, subfolders, source_files):
    subfolders_str = ''.join([f"*   `{d}`\n" for d in subfolders])
    sources_str = ''.join([f"*   `{s}`\n" for s in source_files])
    return f"""
This file has been created automatically. Please adapt it to your needs.

## Content

The content of this DUNE module has been extracted from the module `{module_dir}`.
In particular, the following subfolders of `{module_dir}` have been extracted:
{subfolders_str}
Additionally, all headers in `{module_dir}` that are required to build the
executables from the sources
{sources_str}
have been extracted. You can configure the module just like any other DUNE
module by using `dunecontrol`. For building and running the executables,
please go to the build folders corresponding to the sources listed above.\n
"""


###################################################################
# Installation part of README.md
###################################################################
def info_readme_installation(remoteurl, install_script_name):
    return f"""

## Installation

The easiest way of installation is to use the install script `{install_script_name}`
provided in this repository.
Using `wget`, you can simply install all dependent modules by typing:

```sh
wget {remoteurl}/{install_script_name}.sh
chmod u+x {install_script_name}
./{install_script_name}
```

This will create a sub-folder `DUMUX`, clone all modules into it.

"""


###################################################################
# Infos on how to create the install script manually
###################################################################
def info_make_install(new_module_name):
    return f"""
========================================================================

The extracted module is contained in the subfolder '{new_module_name}'.
You can configure it with
./dune-common/bin/dunecontrol --opts=dumux/cmake.opts —only={new_module_name} all
To create an install script use the "makeinstallscript.py" script in the same folder.
You can run 'python3 makeinstallscript.py --help' for more information
"""


if __name__ == "__main__":

    # set script parameters
    epilog = '''
    -----------------------------------------------------------
    The script has to be called one level above module_dir.
    At least one of the subfolders (FOLDER_1 [FOLDER_2 ...]) has
    to contain a source file *.cc of an executable for which
    you would like to timber a table in dumux-pub.)
    -----------------------------------------------------------

    Example usage:
    ./dumux/bin/extractmodule/extract_as_new_module.py dumux-fracture appl test

    (extracts the subfolders appl and test from the module dumux-fracture)

    '''
    parser = argparse.ArgumentParser(
        prog='extract_as_new_module.py',
        usage='./dumux/bin/extractmodule/extract_as_new_module.py'
              ' module_dir SUBFOLDER_1 [SUBFOLDER_2 ...]',
        description='This script extracts subfolders of a given DUNE module'
                    ' into a new DUNE module.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog)
    parser.add_argument(
        'module_dir',
        help='Dune module from which the subfolder is extracted'
    )
    parser.add_argument(
        'subfolder', nargs='+',
        help='subfolder(s) of module_dir that you want to extract'
    )
    args = vars(parser.parse_args())

    # if module_dir ends with slash(es) remove it/them
    module_dir = args['module_dir'].strip(os.sep)
    module_path = os.path.abspath(module_dir)
    # make unique set of subfolders
    subfolders = list(set(args['subfolder']))

    # check if we are above module_dir
    if not os.path.isdir(module_dir):
        sys.exit("ERROR: You need to run the script"
                 f"one level above the folder {module_dir}.\n"
                 f"Run \"{os.path.basename(__file__)} --help\" for details.")

    # determine all source files in the paths passed as arguments
    source_files = extract_sources_files(module_dir, subfolders)

    # check if sources have been obtained
    if not source_files:
        sys.exit(
            "ERROR: No source files *.cc found in the subfolders: " +
            ", ".join([str(x) for x in subfolders]) + ".\n"
            "Be sure to provide a list of paths as arguments.\n"
            f"Run '{os.path.basename(__file__)} --help' for details."
        )

    # try to find the duneproject script
    dune_project = shutil.which('duneproject', path="dune-common/bin")
    if dune_project is None:
        sys.exit(
            "ERROR: Could not find duneproject.\n"
            "Make sure to have duneproject in dune-common/bin"
        )

    # give explanations
    print(info_explanations(
        module_dir, module_path, subfolders, source_files
    ))
    input("Read the above and press [Enter] to proceed...")

    # run duneproject
    subprocess.call([dune_project])

    # find the created folder
    # as the one with the most recent modification time
    new_module_name = max(
        [d for d in os.listdir() if os.path.isdir(d)],
        key=os.path.getmtime
    )

    # verify it's really a Dune module
    check_module(new_module_name)
    print(
        f"Found new module {new_module_name}\n"
        "Copying source files..."
    )

    # get the base path of the new module
    new_module_path = module_path.replace(module_dir, new_module_name, 1)

    # copy the source tree
    # copy all base folders complete, then delete the unnecessary ones
    base_folders = list(set([
        os.path.relpath(s, module_path).split(os.path.sep)[0] for s in source_files
    ]))
    for b in base_folders:
        path_in_old_module = os.path.join(module_path, b)
        path_in_new_module = os.path.join(new_module_path, b)
        copy_tree(path_in_old_module, path_in_new_module)

    # add base folders in project-level CMakeLists.txt
    with open(os.path.join(new_module_path, "CMakeLists.txt"), "r") as cml:
        section = 0
        content = []
        for line in cml.readlines():
            if line.startswith("add_subdirectory"):
                section = 1 if section in [0, 1] else 2
            elif section == 1:
                section = 2
                for b in base_folders:
                    content.append(f"add_subdirectory({b})")
            content.append(line)

    with open(os.path.join(new_module_path, "CMakeLists.txt"), "w") as cml:
        for line in content:
            cml.write(line)

    # go through source tree and remove unnecessary directories
    new_source_files = [
        s.replace(module_dir, new_module_name, 1) for s in source_files
    ]
    for b in base_folders:
        def matching_path(path, list_of_paths):
            for p in list_of_paths:
                if path in p:
                    return True
            return False

        for path, dirs, files in os.walk(os.path.join(new_module_path, b)):
            for d in dirs:
                dir_path = os.path.join(path, d)
                keep = matching_path(dir_path, new_source_files)
                if not keep:
                    rel_dir_path = os.path.relpath(dir_path, new_module_path)
                    answer_is_yes = query_yes_no(
                        f"{rel_dir_path} does not contain source files."
                        " Copy to new module?",
                        default="no"
                    )
                    if not answer_is_yes:
                        # remove copy of directory
                        shutil.rmtree(dir_path)
                        # remove entry from CMakeLists.txt
                        cml_path = os.path.join(module_path, "CMakeLists.txt")
                        with open(cml_path, "r+") as cml:
                            content = cml.read()
                            cml.seek(0)
                            content = content.replace(
                                "add_subdirectory({d})\n", ""
                            )
                            cml.write(content)
                            cml.truncate()

                        # make os.walk know about removed folders
                        dirs.remove(d)

    # search for all header (in parallel)
    with mp.Pool() as p:
        headers = itertools.chain.from_iterable(p.map(
            partial(search_headers, module_path=module_path),
            source_files
        ))

    # make unique
    headers = list(set(headers))

    # copy headers to the new module
    for header in headers:
        header_dir = os.path.dirname(os.path.realpath(header))
        path_in_new_module = header_dir.replace(module_dir, new_module_name, 1)
        os.makedirs(path_in_new_module, exist_ok=True)
        shutil.copy(header, path_in_new_module)

    # copy .gitignore from dumux to the new module
    dumux_gitignore_file = "dumux/.gitignore"
    shutil.copy(dumux_gitignore_file, new_module_path)

    # delete unnecessary directories to keep the extracted module clean
    if "dune" not in subfolders:
        shutil.rmtree(os.path.join(new_module_path, 'dune'))
    if "src" not in subfolders:
        shutil.rmtree(os.path.join(new_module_path, 'src'))
    with open(os.path.join(new_module_path, "CMakeLists.txt"), "r+") as cml:
        content = cml.read()
        cml.seek(0)
        if "dune" not in subfolders:
            content = content.replace("add_subdirectory(dune)\n", "")
        if "src" not in subfolders:
            content = content.replace("add_subdirectory(src)\n", "")
        cml.write(content)
        cml.truncate()

    # create README file
    os.remove(os.path.join(new_module_path, 'README'))
    readme_path = os.path.join(new_module_path, "README.md")
    with open(readme_path, "w") as readme_file:
        readme_file.write(
            info_readme_main(module_dir, subfolders, source_files)
        )

    # ask user if to write version information into README.md
    if query_yes_no("Write detailed version information"
                    " (folder/branch/commits/dates) into README.md?"):
        print("Looking for the dune modules in path: "
              + os.path.abspath(".")
              + "...")

        versions = getPersistentVersions(
            [dep['folder'] for dep in getDependencies(module_path)], True
        )

        # append version information
        with open(readme_path, "a") as readme_file:
            readme_file.write(
                "\n## Version Information\n\n" +
                versionTable(versions)
            )

    # if there is a remote repository available
    # we can directly push the source code and
    # also create a "one-click" installation script
    install_script_name = 'install_' + new_module_name + '.sh'
    if query_yes_no("Do you have an empty remote repository to push the code to (recommended)?"):
        remoteurl = get_remote_url(new_module_path)

        run_from_mod = callFromPath(new_module_path)(runCommand)
        try:
            run_from_mod('git init')
            run_from_mod('git add .')
            run_from_mod('git commit -m "Initial commit"')
            run_from_mod('git remote add origin {}'.format(remoteurl))
            run_from_mod('git push -u origin master')

            # create an installation script
            makeInstallScript(new_module_path, ignoreUntracked=True)
            shutil.move(install_script_name,
                        os.path.join(new_module_path, install_script_name))

            # append install information into readme
            with open(readme_path, "a") as readme_file:
                readme_file.write(info_readme_installation(remoteurl, install_script_name))

            # add to version control
            run_from_mod('git add .')
            run_from_mod('git commit -m "Create Install script"')
            run_from_mod('git push -u origin master')

        except Exception:
            sys.exit(
                f"Automatically generate install script {install_script_name} failed."
                "\nTo create the script, use the 'makeinstallscript.py' script in the same folder."
                "\nRun 'python3 makeinstallscript.py --help' for more detailed information"
            )

    # output guidance for users to create install script manually
    else:
        print(info_make_install(new_module_name))
