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
import threading
import fnmatch
from util import getPersistentVersions
from util import printVersionTable
from makeinstallscript import makeInstallScript
try:
    path = os.path.split(os.path.abspath(__file__))[0]
    sys.path.append(os.path.join(path, '../bin/util'))
    from getmoduleinfo import getDependencies
    from common import callFromPath, runCommand
except Exception:
    sys.exit('Could not import common modul or getModuleInfo')

# set script parameters
epilog = '''
-----------------------------------------------------------
The script has to be called one level above module_dir.
At least one of the subfolders (FOLDER_1 [FOLDER_2 ...]) has
to contain a source file *.cc of an executable for which
you would like to timber a table in dumux-pub.)
'''
parser = argparse.ArgumentParser(
    prog='extractmodulepart',
    usage="./extractmodulepart module_dir SUBFOLDER_1 [SUBFOLDER_2 ...]",
    description='This script extracts a subfolder of a DUNE module',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=epilog)
parser.add_argument('module_dir',
                    help='Dune module from which the subfolder is extracted')
parser.add_argument('subfolder',
                    nargs='+',
                    help='subfolder(s) of module_dir that you want to extract')
args = vars(parser.parse_args())


# function to search matching header files
def search_headers(c_file):
    with open(c_file, 'r') as f:
        content = f.read()

    header_in_bracket = re.findall(r'#include\s+<(.+?)>', content)
    header_in_quotation = re.findall(r'#include\s+"(.+?)"', content)

    # search for includes relative to the module path
    for header in header_in_bracket + header_in_quotation:
        header_with_path = os.path.join(module_full_path, header)
        if not os.path.exists(header_with_path):
            continue
        if header_with_path not in all_headers:
            all_headers.append(header_with_path)
            thread_ = threading.Thread(target=search_headers,
                                       args=(header_with_path, ))
            thread_.start()

    # search for includes relative to the path of the including file
    # only allow quoted includes for this
    for header in header_in_quotation:
        if header == "config.h":
            continue
        header_dir_name = os.path.dirname(c_file)
        header_with_path = os.path.join(header_dir_name, header)
        if not os.path.exists(header_with_path):
            continue
        if header_with_path not in all_headers:
            all_headers.append(header_with_path)
            thread_ = threading.Thread(target=search_headers,
                                       args=(header_with_path, ))
            thread_.start()


# functions to find the file with specific name in a direcotry
def find_files_with_name(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result


# functions to find the file with specific pattern in a direcotry
def find_files_with_pattern(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


# functions to configure CMakeLists.txt files
def generate_new_content(lines, pattern, replace_str):
    index = lines.find(pattern)
    if index != -1:
        content = lines[0: index] + "\n"
        if replace_str != "":
            content += replace_str + "\n"
        flag = True
        while flag:
            index = lines.find(")", index)
            if index != -1:
                lines = lines[index + 1:]
            else:
                break

            index = lines.find(pattern)
            if index != -1:
                content += lines[0: index] + "\n"
            else:
                content += lines + "\n"
                flag = False
    else:
        if replace_str == "":
            content = lines
        else:
            if pattern == "add_subdirectory(":
                content = replace_str + "\n" + lines
            else:
                content = lines + "\n" + replace_str
    return content


def generate_subdirectory_content(dirs):
    content = ""
    for name in dirs:
        content += "add_subdirectory(" + name + ")" + "\n"
    return content


def generate_install_content(header_files, destination):
    if len(header_files) == 0:
        return ""
    content = "install(FILES" + "\n"
    for header_file in header_files:
        content += "    " + header_file + "\n"
    content += ("DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}"
                + destination + ")" + "\n")
    return content


def generate_content(cmake_list_txt_file, dirs, header_files, destination):
    subdirectory_content = generate_subdirectory_content(dirs)
    install_content = generate_install_content(header_files, destination)

    if os.path.exists(cmake_list_txt_file):
        with open(cmake_list_txt_file, "r", encoding="utf-8") as f:
            content = "".join(f.readlines()).strip()
            f.close()
        content = generate_new_content(content, "add_subdirectory(",
                                       subdirectory_content)
        return generate_new_content(content, "install(FILE", install_content)
    else:
        if subdirectory_content == "":
            return install_content
        else:
            return subdirectory_content + "\n" + install_content


def generate_cmake_lists_txt(cmake_list_txt_file, dirs,
                             header_files, destination):
    content = generate_content(cmake_list_txt_file, dirs,
                               header_files, destination)

    with open(cmake_list_txt_file, "w", encoding="utf-8") as f:
        if content != "":
            f.write(content)


def check_dir(root_dir):
    if not os.path.exists(root_dir):
        print("root path" + str(root_dir) + "is not exist!")
        return False
    if not os.path.isdir(root_dir):
        print("root path" + str(root_dir) + "is not dir!")
        return False
    return True


def check_str(root_dir):
    if root_dir is None or root_dir.strip() == "":
        return None
    root_dir = root_dir.strip()
    if root_dir.endswith(os.sep):
        root_dir = root_dir[:-1]
    return root_dir


def generate_cmake_lists_txt_file(root_dir):
    root_dir = check_str(root_dir)
    if root_dir is None:
        return
    if not check_dir(root_dir):
        return
    drop_len = len(root_dir)
    for parent, dirs, files in os.walk(root_dir):
        destination = parent[drop_len:]
        header_files = []
        cmake_list_txt_file = os.path.join(parent, "CMakeLists.txt")
        for name in files:
            for suffix in [".h", ".hh"]:
                if name.endswith(suffix):
                    header_files.append(name)
        generate_cmake_lists_txt(cmake_list_txt_file, dirs,
                                 header_files, destination)


# function asking user to answer yes or no question
def query_yes_no(question, default="yes"):
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [y/N] "
    elif default == "no":
        prompt = " [Y/n] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)
    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == "":
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


# funtion to get correct git remote url and make sure it is empty
def get_remote_url(repo_path):
    run_from_mod = callFromPath(repo_path)(runCommand)
    while True:
        if query_yes_no("Do you own a subfolder in dumux-pub?"):
            nameyearx = input("Type below the information of"
                              " AuthorLastNameYearx (e.g. luggi2020b)\n")
            remoteurl = 'https://git.iws.uni-stuttgart.de'\
                        '/dumux-pub/{}.git'.format(nameyearx)
        else:
            remoteurl = input("Provide URL of your remote repository:\n")
        check_remote_repo = run_from_mod('git ls-remote {}'.format(remoteurl))
        if (check_remote_repo is None):
            sys.stdout.write("\nERROR: Please re-enter correct information.\n")
        elif (check_remote_repo == ''):
            return remoteurl
        else:
            sys.stdout.write("ERROR: The remote reposity is not empty!.\n")


if __name__ == "__main__":

    module_dir = args['module_dir']
    subfolders = args['subfolder']

    # if module_dir contains a slash as last character, delete it
    if module_dir.endswith('/'):
        module_dir = module_dir[:-1]

    # check if we are above module_dir
    if not (os.path.isdir(module_dir)):
        print("ERROR: you need to run the script one level above the folder "
              + module_dir + ".")
        print("Run \""+os.path.basename(__file__)+" --help\" for details.")
        exit(1)

    # determine all source files in the paths passed as arguments
    script_path = os.getcwd()
    os.chdir(module_dir)
    module_full_path = os.getcwd()
    all_sources = []
    all_sources_with_path = []
    all_directories = []
    for dir_path in subfolders:
        if dir_path.startswith(module_dir):
            # the "+1" is required to also remove the "/"
            stripped_path = dir_path[len(module_dir)+1:]
        else:
            stripped_path = dir_path
        full_path = os.path.join(module_full_path, stripped_path)
        directories = " " + stripped_path
        all_directories.append(stripped_path)
        os.chdir(full_path)
        for file in find_files_with_pattern("*.cc", full_path):
            print(file)
            sources = os.path.relpath(file, module_full_path)
            print(sources)
            all_sources.append(sources)
            sourceswithpath = os.path.join(module_full_path, sources)
            all_sources_with_path.append(sourceswithpath)
    os.chdir(module_full_path)
    os.chdir("..")  # back to the script folder

    # check if sources have been obtained
    if (all_sources == []):
        print("ERROR: no source files *.cc found in the directories "
              + " ".join([str(x) for x in subfolders]) + ".")
        print("Be sure to provide a list of paths as arguments.")
        print("Run \""+os.path.basename(__file__)+" --help\" for details.")
        exit(1)

    # try to find the duneproject script
    dune_project = shutil.which('duneproject', path="dune-common/bin")
    if (dune_project is None):
        print("ERROR: Could not find duneproject.")
        print("Be sure to have duneproject in dune-common/bin")
        exit(1)
    else:
        print(dune_project)

    # give explanations
    print("\n""This script will\n"
          "- extract the following sub-folders of " + module_dir + ":\n")
    for dir_path in all_directories:
        print("  " + dir_path + ",")
    print("\n""  and all headers in " + module_dir
          + " that are required to build the\n"
          "  executables from the sources\n")
    for source in all_sources:
        print("  " + source + ",")
    print("""
This script copies extracted files into a freshly created DUNE module
retaining the directory structure
and updates/creates required files like CMakeLists.txt.
duneproject will run now.
The new module should NOT depend on the module in {0}.\n\n"""
          .format(module_dir))

    input("Read the above and press [Enter] to proceed...")

    # run duneproject
    old_ls = os.listdir()
    subprocess.call([dune_project])
    new_ls = os.listdir()

    # determine the new module/directory name
    module_name = (set(new_ls) - set(old_ls)).pop()
    if (module_name == ""):
        print("ERROR: could not find new module. Aborting.")
        exit(1)
    else:
        print()
        print(os.path.basename(__file__) + ": Found new module " + module_name)
    print("Determining required headers...")
    os.chdir(module_name)
    module_path = os.getcwd()

    # extract headers, provide some output and copy everything to new module
    all_headers = []
    print("The following header files are extracted: ")
    for source in all_sources_with_path:
        search_headers(source)
        for header in all_headers:
            print(header)
            dir_path = os.path.dirname(os.path.realpath(header))\
                .replace(module_dir, module_name, 1)
            os.makedirs(dir_path, exist_ok=True)
            shutil.copy(header, dir_path)
        source_dir = os.path.dirname(source)
        source_path = source_dir.replace(module_dir, module_name, 1)
        copy_tree(source_dir, source_path)

    # copy .gitignore from dumux to the new module
    dumux_gitignore_file = os.path.join(script_path, "dumux/.gitignore")
    shutil.copy(dumux_gitignore_file, module_path)

    # delete unnecessary directories to keep the extracted module clean
    shutil.rmtree('dune')
    shutil.rmtree('src')

    # set CMakeLists.txt for each directory
    generate_cmake_lists_txt_file(module_path)
    print("The required header files are extracted"
          " and CMakeLists are configured.")
    print("==================================================================")

    # create README file
    os.remove("README")
    orig_stdout = sys.stdout
    readme_file = open("README.md", "w+")

    readme_dir_str_lists = []
    for dir_path in all_directories:
        readme_dir_str_lists.append("*   `" + dir_path + "`,\n")
    readme_source_str_lists = []
    for source in all_sources:
        readme_source_str_lists.append("*   `" + source + "`,\n")
    readme_str_args = [module_dir, ''.join(readme_dir_str_lists),
                       ''.join(readme_source_str_lists)]

    readme_file.write("""This file has been created automatically. Please adapt it to your needs. \n
===============================\n
## Content\n
"The content of this DUNE module was extracted from the module `{0}`."
In particular, the following subfolders of `{0}` have been extracted:
{1}
Additionally, all headers in `{0}` that are required to build the
executables from the sources
{2}
have been extracted. You can configure the module just like any other DUNE
module by using `dunecontrol`. For building and running the executables,
please go to the build folders corresponding to the sources listed above.\n
""".format(*readme_str_args))

    version_script_path = os.path.join(script_path,
                                       "dumux/bin/util/getusedversions.py")
    install_script_path = os.path.join(script_path,
                                       "dumux/bin/util/makeinstallscript.py")
    os.chdir(script_path)
    # ask user if to write version information into README
    if query_yes_no("\nWrite detailed version information"
                    " (folder/branch/commits/dates) into README.md?"):
        print("Looking for the dune modules in path :" + str(script_path)
              + "...")
        readme_file.write("===============================\n"
                          + "## Version Information\n")
        versions = getPersistentVersions(
            [dep['folder'] for dep in getDependencies(module_full_path)], True)
        print("Writing version information into README.md ...")
        sys.stdout = readme_file
        printVersionTable(versions)
        sys.stdout = orig_stdout

    print("Automatic generation of README.md file is complete.")
    print("==================================================================")

    # ask user if to generate an install script
    install_script_name = 'install_' + module_name + '.sh'
    if query_yes_no("\nGenerate install script " + install_script_name
                    + " in your new module " + module_name + "?"):
        remoteurl = get_remote_url(module_path)
        print("Write instructions for install script into README.md file...")
        readme_file.write("""===============================
## Installation
   The easiest way of installation is to use the install script `{0}`
   provided in this repository.
   Using `wget`, you can simply install all dependent modules by typing:

   ```sh
   wget {1}/{0}.sh
   chmod u+x {0}
   ./{0}
   ```

   This will create a sub-folder `DUMUX`, clone all modules into it.

""".format(*[install_script_name, remoteurl]))
        readme_file.close()
        run_from_mod = callFromPath(module_path)(runCommand)
        try:
            run_from_mod('git init')
            run_from_mod('git add .')
            run_from_mod('git commit -m "Initial commit"')
            run_from_mod('git remote add origin {}'.format(remoteurl))
            run_from_mod('git push -u origin master')

            makeInstallScript(os.path.join(script_path, module_name),
                              ignoreUntracked=True)
            shutil.move(install_script_name,
                        os.path.join(module_name, install_script_name))
            run_from_mod('git add .')
            run_from_mod('git commit -m "Create Install script"')
            run_from_mod('git push -u origin master')

            print("\n" + "*"*80 + "\n"
                  """Congratulations, the install script {}.sh has been succesfully generated,
                  You can make further changes to suit your needs,
                  commit and push to the remote repository"""
                  .format(module_name) + "\n" + "*"*80)
        except Exception:
            sys.exit("""Automatically generate install script {}.sh failed.
To create the script, use the "makeinstallscript.py" script in the same folder,
run python3 makeinstallscript.py --help for more detailed information"""
                     .format(module_name))

    # in case user want to do generate install.sh manuelly
    else:
        # output guidence for users
        print("\n" + "*"*80 + "\n" +
              """The extracted module is contained in the subfolder \"{0}\".
              You can configure it using \"dunecontrol ... —only= {0} all\".
              Make sure you have an empty remote repository.
              To create an install script {0}.sh,
              use the "makeinstallscript.py" script in the same folder,
              run python3 makeinstallscript.py --help for more information"""
              .format(module_name) + "\n" + "*"*80)
