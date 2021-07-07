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
from pathlib import Path
from functools import partial
from util import getPersistentVersions
from util import versionTable
from makeinstallscript import makeInstallScript
from makeinstallscript import python_or_bash
try:
    path = os.path.split(os.path.abspath(__file__))[0]
    sys.path.append(os.path.join(path, '../bin/util'))
    from getmoduleinfo import getDependencies
    from common import callFromPath, runCommand
except Exception:
    sys.exit('Could not import common modul or getModuleInfo')
import logging


# return the list of included headers including the header itself
def add_headers_recursively(header_path, current_headers, module_path):
    if os.path.exists(header_path):
        if header_path not in current_headers:
            current_headers.append(header_path)
            search_headers(header_path, current_headers, module_path)


# function to search matching header files
def search_headers(source_file, headers, module_path):
    with open(source_file, 'r') as f:
        content = f.read()
        header_in_bracket = re.findall(r'#include\s+<(.+?)>', content)
        header_in_quotation = re.findall(r'#include\s+"(.+?)"', content)

    # search for includes relative to the module path
    for header in header_in_bracket:
        header_path = os.path.join(module_path, header)
        add_headers_recursively(header_path, headers, module_path)

    # search for includes relative to the path of the including file
    # only allow quoted includes for this
    for header in header_in_quotation:
        if header == "config.h":
            continue
        header_path = os.path.join(os.path.dirname(source_file), header)
        add_headers_recursively(header_path, headers, module_path)

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
            logging.error("Remote repository is not empty!")
            sys.stdout.write("ERROR: The remote reposity is not empty!.\n")


def path_check(basedir, subdir):
    # check if basedir contained in the script path
    if not os.path.isdir(basedir):
        logging.error(f"No {basedir} found in your path where the script is running!")
        sys.exit("ERROR: You need to run the script"
                 f"one level above the folder {basedir}.\n"
                 f"Run \"{os.path.basename(__file__)} --help\" for details.")

    # check whether subdir is a sub-directory of basedir
    base = Path(basedir)
    for subfolder in subdir:
        child = Path(os.path.join(base, subfolder))
        if base not in child.parents or not os.path.isdir(child):
            logging.error(f"Subfolder '{subfolder}' is not a subfolder of '{basedir}'")
            raise NameError(f"Subfolder '{subfolder}' is not a subfolder of '{basedir}'")


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
        logging.error("The module you created is not dune-module!")
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
that are required to build the exectutables from the sources:
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
def info_readme_installation(remoteurl, install_script_name, new_module_name, language):
    return f"""

## Installation

The installation procedure is done as follows :
Create a root folder, e.g. `DUMUX`, enter the previously created folder,
clone the remote and use the install script `{install_script_name}`
provided in this repository to install all dependent modules.

```sh
mkdir DUMUX
cd DUMUX
git clone {remoteurl}
{"./" if language == "bash" else "python3 "}{new_module_name}/{install_script_name}
```

This will clone all modules into the directory `DUMUX`,
configure your module with `dunecontrol` and build tests.

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
    # generate a log file to control all the process
    logging.basicConfig(filename="extractmodulepart.log", filemode='w', level=logging.DEBUG)
    logging.info("The main function is going to be executed, preparing to take off.")
    logging.debug("Passing arguments provided by user...")
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
    logging.debug("Parameters from terminal are passed.")

    # if module_dir ends with slash(es) remove it/them
    module_dir = args['module_dir'].strip(os.sep)
    module_path = os.path.abspath(module_dir)
    # make unique set of subfolders
    subfolders = list(set(args['subfolder']))

    # check paths to prevenet possible errors
    logging.debug("Checking if base module and subfolders exist...")
    path_check(module_dir, subfolders)
    logging.debug("Path check is done.")

    # determine all source files in the paths passed as arguments
    source_files = extract_sources_files(module_dir, subfolders)

    # check if sources have been obtained
    logging.debug("Checking if sources have been obtained...")
    if not source_files:
        logging.error("No sources found in the subfolders you provide!")
        sys.exit(
            "ERROR: No source files *.cc found in the subfolders: " +
            ", ".join([str(x) for x in subfolders]) + ".\n"
            "Be sure to provide a list of paths as arguments.\n"
            f"Run '{os.path.basename(__file__)} --help' for details."
        )
    logging.debug("Sources are found in subfolders.")

    # try to find the duneproject script
    dune_project = shutil.which('duneproject', path="dune-common/bin")
    if dune_project is None:
        logging.error("No duneproject found in dune-common/bin!")
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
    try:
        logging.info("Calling dune-common/bin/duneproject to create the new extracted module...")
        subprocess.call([dune_project])
        logging.info("--The new module is created sucessfully.")
    except Exception:
        logging.error("Failed to generate new module with duneproject!")

    # find the created folder
    # as the one with the most recent modification time
    logging.debug("Getting the name of the new module...")
    new_module_name = max(
        [d for d in os.listdir() if os.path.isdir(d)],
        key=os.path.getmtime
    )
    logging.debug(f"Name of the new module is detected as {new_module_name}.")

    # verify it's really a Dune module
    logging.debug("Checking if the new module is dune module...")
    check_module(new_module_name)
    logging.debug("The new module is checked as dune module.")
    print(
        f"Found new module {new_module_name}\n"
        "Copying source files..."
    )

    logging.info("Extracting required headers and copy useful files to the new module...")
    # get the base path of the new module
    logging.debug("Specifying paths of new module...")
    new_module_path = new_module_name.join(module_path.rsplit(module_dir, 1))
    logging.debug(f"The path for new module is {new_module_path}.")

    # copy the source tree
    # copy all base folders complete, then delete the unnecessary ones
    logging.debug("Copying old directories cotaining source files into new module...")
    for b in subfolders:
        path_in_old_module = os.path.join(module_path, b)
        path_in_new_module = os.path.join(new_module_path, b)
        copy_tree(path_in_old_module, path_in_new_module)
    logging.debug("Directoies containing source files are copied.")

    # add base folders in project-level CMakeLists.txt
    logging.debug("Adding base folders to CMakeList files in new module...")
    with open(os.path.join(new_module_path, "CMakeLists.txt"), "r") as cml:
        section = 0
        content = []
        for line in cml.readlines():
            if line.startswith("add_subdirectory"):
                section = 1 if section in [0, 1] else 2
            elif section == 1:
                section = 2
                for b in subfolders:
                    content.append(f"add_subdirectory({b})")
            content.append(line)

    with open(os.path.join(new_module_path, "CMakeLists.txt"), "w") as cml:
        for line in content:
            cml.write(line)
    logging.debug("CMakelist at top level in new module is confiugred.")

    # search for all header (in parallel)
    logging.debug("Searching for all header files...")
    with mp.Pool() as p:
        headers = itertools.chain.from_iterable(p.map(
            partial(search_headers, module_path=module_path, headers=[]),
            source_files
        ))
    logging.debug("Head files are found.")

    # make unique
    logging.debug("Making header files unqiue (removing duplicates)...")
    headers = list(set(headers))
    logging.debug("Duplicates are removed.")

    # copy headers to the new module
    logging.debug("Copying headers to the new module...")
    header_dirs = []
    for header in headers:
        header_dir = os.path.dirname(os.path.realpath(header))
        path_in_new_module = header_dir.replace(module_dir, new_module_name, 1)
        header_dirs.append(path_in_new_module)
        os.makedirs(path_in_new_module, exist_ok=True)
        shutil.copy(header, path_in_new_module)
    logging.debug("Headers are copied to the new module.")

    # copy .gitignore from dumux to the new module
    dumux_gitignore_file = "dumux/.gitignore"
    logging.debug("Copying .git from dumux into new module...")
    shutil.copy(dumux_gitignore_file, new_module_path)
    logging.debug(".gitignore is copied.")
    logging.info("--Requried headers are extracted and useful files are copied into the new module succesfully.")

    # delete unnecessary directories to keep the extracted module clean
    logging.info("Cleaning up the new module...")
    # go through source tree and remove unnecessary directories
    logging.debug("Handling directories where no source is obtained...")
    new_source_files = [
        s.replace(module_dir, new_module_name, 1) for s in source_files
    ]
    for b in subfolders:
        def matching_path(path, list_of_paths, header_dirs):
            for p in list_of_paths:
                if path in p or path in header_dirs:
                    return True
            return False

        no_source_folder = []
        for path, dirs, files in os.walk(os.path.join(new_module_path, b)):
            for d in dirs:
                dir_path = os.path.join(path, d)
                keep = matching_path(dir_path, new_source_files, header_dirs)
                if not keep:
                    rel_dir_path = os.path.relpath(dir_path, new_module_path)
                    no_source_folder.append(rel_dir_path)

        no_source_lists = ('\n'.join(dir for dir in no_source_folder))
        yes_to_all = query_yes_no(
            "Could not automatically determine if following directories contain data essential for the extracted applications:\n{0}.\n"
            "Do you want to copy all of them in the new module {1}?\n"
            "(By choosing no you need to consider each folder seperately)"
            .format(*[no_source_lists, new_module_name]),
            default="yes"
        )

        if not yes_to_all:
            for rel_dir_path in no_source_folder:
                answer_is_yes = query_yes_no(
                    f"Could not automatically determine if {rel_dir_path} contains data essential for the extracted applications.\n"
                    f"Do you want to copy the files in {rel_dir_path} to the extracted module {new_module_name}",
                    default="yes"
                )
                if not answer_is_yes:
                    # remove copy of directory
                    shutil.rmtree(os.path.join(new_module_path, rel_dir_path))
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
    logging.debug("The folders containg no sources are handled.")

    logging.debug("Deleting dune/src in the extracted module...")
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
    logging.debug("Dune/src are deleted.")
    logging.info("--The new module is cleaned up successfully.")

    # create README file
    logging.info("Creating README file and generating an install script in new module...")
    os.remove(os.path.join(new_module_path, 'README'))
    readme_path = os.path.join(new_module_path, "README.md")
    with open(readme_path, "w") as readme_file:
        readme_file.write(
            info_readme_main(module_dir, subfolders, source_files)
        )
    logging.debug("README file in new module is created.")

    # ask user if to write version information into README.md
    logging.debug("Ask user if to write version information into README file.")
    if query_yes_no("Write detailed version information"
                    " (folder/branch/commits/dates) into README.md?\n"):
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
        logging.debug("The version information is written into README file.")

    # if there is a remote repository available
    # we can directly push the source code and
    # also create a "one-click" installation script
    logging.debug("Generating an install script...")
    logging.debug("Ask user to choose python or bash to generate the install script.")
    language = python_or_bash()
    install_script_name = 'install_' + new_module_name + '.%s'%("sh" if language == "bash" else "py")
    try:
        makeInstallScript(new_module_path, ignoreUntracked=True, skipFolders=new_module_name,
                          suppressHints=True, topFolderName=None, language=language)
        logging.debug("The install script is generated.")
    except Exception:
        logging.error(f"Failed to generate script {install_script_name} by calling external functino makeInstallScript!")
        sys.exit(
            f"Automatically generate install script {install_script_name} failed."
            "\nTo create the script, use the 'makeinstallscript.py' script in the same folder."
            "\nRun 'python3 makeinstallscript.py --help' for more detailed information"
        )

    shutil.move(install_script_name,
                os.path.join(new_module_path, install_script_name))
    print(info_make_install(new_module_name))
    logging.info("--README file and install script are generated succesfully.")
    run_from_mod = callFromPath(new_module_path)(runCommand)

    logging.info("Commiting the new module and pushing to remote if url is provided...")
    if query_yes_no("Do you have an empty remote repository to push the code to (recommended)?"):
        logging.debug("Trying to get the remote URL.")
        remoteurl = get_remote_url(new_module_path)
        logging.debug(f"The remote URL proviede by the user is {remoteurl}")

        # append install information into readme
        logging.debug("Writing README file of the new module...")
        with open(readme_path, "a") as readme_file:
            readme_file.write(info_readme_installation(remoteurl, install_script_name, new_module_name, language))
        logging.debug("The README file is written.")
        logging.debug("Initializing git repository...")
        run_from_mod('git init')
        logging.debug("Git repository is initialized.")
        logging.debug("Adding all the files into git staging area...")
        run_from_mod('git add .')
        logging.debug("All files are added into staging area.")
        logging.debug("Commiting all files in staging area as inital commit...")
        run_from_mod('git commit -m "Initial commit"')
        logging.debug("All files in staging area are commited.")
        logging.debug(f"Adding remote {remoteurl} to the origin...")
        run_from_mod('git remote add origin {}'.format(remoteurl))
        logging.debug(f"{remoteurl} added to git remote origin.")
        logging.debug(f"Pushing all the changes into remote git repository {remoteurl}")
        run_from_mod('git push -u origin master')
        logging.debug("Commits are pushed to git remote.")

    # output guidance for users to create install script manually
    else:
        remoteurl = "{$remoteurl$} (needs to be manuelly adapted later)"
        # append install information into readme
        logging.debug("Writing README file of the new module...")
        with open(readme_path, "a") as readme_file:
            readme_file.write(info_readme_installation(remoteurl, install_script_name, new_module_name, language))
        logging.debug("The README file is written.")
        logging.debug("Initializing git repository...")
        run_from_mod('git init')
        logging.debug("Git repository is initialized.")
        logging.debug("Adding all the files into git staging area...")
        run_from_mod('git add .')
        logging.debug("All files are added into staging area.")
        logging.debug("Commiting all files in staging area as inital commit...")
        run_from_mod('git commit -m "Initial commit"')
        logging.debug("All files in staging area are commited.")
        print("\nPlease remember to replace placeholder $remoteurl$ in installation part\n"
              "of the README file with the correct URL of your remote git repository.")
    logging.info("--Changes are commited (and pushed if remote is known) successfully.")
    logging.info("Congratulations, everything is fine and landing is smooth!")
