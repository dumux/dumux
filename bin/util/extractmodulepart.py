#!/usr/bin/env python3
import sys, glob, os, subprocess
import argparse
import shutil
from distutils.dir_util import copy_tree
import os.path

"""
This is a python script for extracting modules.
It was originally written to solve the problem that the bash shell script could only run in a specific environment.
So the whole design philosophy is to maximize the portability, i.e. write in a cross-platform style.
"""

def show_helpmesseage():
    print("\n""USAGE: "+os.path.basename(__file__)+" module_dir FOLDER_1 [FOLDER_2 ...]\n\n"
          "module_dir is the folder containing the DUNE module from which you\n"
          "want to extract. The script has to be called one level above it.\n\n"
          "The FOLDERs need to indicate subfolders of module_dir. At least one\n"
          "of them has to contain a source file *.cc of an executable for which\n"
          "you would like to timber a table in dumux-pub.")
    exit(1)

# check if help is needed
if (len(sys.argv) < 2):
    show_helpmesseage()
if (str(sys.argv[1]) == "--help" or
    str(sys.argv[1]) == "-help" or
    str(sys.argv[1]) == "help"):
    show_helpmesseage()
module_dir = str(sys.argv[1])

# if module_dir contains a slash as last character, delete it
if module_dir.endswith('/'):
   module_dir = module_dir[:-1]

# check if we are above module_dir
if not (os.path.isdir(module_dir)):
    print("ERROR: you need to run the script one level above the folder "+module_dir+".")
    print("Run \""+os.path.basename(__file__)+" --help\" for details.")
    exit(1)

# determine all source files in the paths passed as arguments
os.chdir(module_dir)
module_full_path=os.getcwd()
all_sources=[]
all_directories=[]
for dir_path in sys.argv[2:]:
    stripped_path = dir_path.removeprefix(module_dir)
    directories = " " + stripped_path
    all_directories.append(stripped_path)
    os.chdir( os.path.join(module_full_path + stripped_path))
    for file in glob.glob("*.cc"):
        sources = os.path.join(stripped_path, file)
        all_sources.append(sources)
os.chdir(module_full_path)
os.chdir("..")

# check if sources have been obtained
for source in all_sources:
    contracted = str(os.popen("echo \"" + source + "\" | tr -d \" tnr\"").read().rstrip("\n"))
    if ( contracted == "" ):
        print("ERROR: no source files *.cc found in the directories "+sys.argv[2]+".")
        print("Be sure to provide a list of paths as arguments to this script.")
        print("Run \""+os.path.basename(__file__)+" --help\" for details.")
        exit(1)

# try to find the duneproject script
dune_project = shutil.which('duneproject', path = "dune-common/bin")
if (dune_project == None):
    print("ERROR: Could not find duneproject.")
    print("Be sure to either have duneproject in your search path")
    print("or to run this script from a directory that contains duneproject.")
    exit(1)
else:
    print(dune_project)

# give explanations
print("\n""This script will\n"
      "- extract the following sub-folders of "+ module_dir +":\n\n")
for dir_path in all_directories:
    print("  "+ dir_path + ",")
print("\n""  and all headers in "+ module_dir + " that are required to build the\n"
      "  executables from the sources\n\n")
for source in all_sources:
    print("  "+ source +",")
print("\n- copy the extracted files into a freshly created DUNE module, retaining the\n"
      "  directory structure,\n\n"
      "- update/create all required files like CMakeLists.txt,\n\n"
      "- store the versions of all used Dune module\n\n"
      "- and extract their modifications as patches.\n\n"
      "Thus, you receive a fully-working DUNE module containing the subset of\n"
       + module_dir + " that is required to run your application.\n"
      "duneproject will be run now. The new module should NOT depend on the\n"
      "module in "+ module_dir + ".\n\n"
      "Read the above and press [Enter] to proceed...",end="")

# run duneproject
old_ls = os.listdir()
subprocess.call([dune_project], shell=True)
new_ls  = os.listdir()

# determine the new module/directory name
module_name = (set(new_ls) - set(old_ls)).pop()
if ( module_name  == "" ):
    print("ERROR: could not find new module. Aborting.")
    exit(1)
else:
    print()
    print( os.path.basename(__file__) + ": Found new module " + module_name)
print("Determining required headers...",end="")
os.chdir(module_name)
module_path=os.getcwd()

# extract all headers
import re
import threading

all_headers = []
all_headers_tmp = []
def search_headers(c_file):
    global all_headers
    f = open(c_file,'r')
    content = f.read()
    f.close()
    header_in_bracket = re.findall(r'(?<=#include <).+?(?=>)',content)
    header_in_quotation = re.findall(r'(?<=#include ").+?(?=")',content)
    for header in header_in_bracket:
        if header.startswith("dumux"):
            header_with_path = os.path.join(module_full_path, header)
        else:
            continue
        if header_with_path not in all_headers:
            all_headers.append(header_with_path)
            thread_ = threading.Thread(target = search_headers, args = (header_with_path, ))
            thread_.start()
    for header in header_in_quotation:
        header_dir_name = os.path.dirname(c_file)
        header_with_path = os.path.join(header_dir_name, header)
        if header_with_path not in all_headers:
            all_headers.append(header_with_path)
            thread_ = threading.Thread(target = search_headers, args = (header_with_path, ))
            thread_.start()

for source in all_sources:
    source_abspath = module_full_path+source
    search_headers(source_abspath)
    for header in all_headers:
        print(header)
        dir_path = os.path.dirname(os.path.realpath(header)).replace("dumux",module_name,1)
        os.makedirs(dir_path, exist_ok=True)
        shutil.copy(header,dir_path)
    source_dir = source_abspath.removesuffix("main.cc")
    source_path = source_dir.replace("dumux",module_name,1)
    copy_tree(source_dir, source_path)

# delete all architecture-dependent files and unneeded directories
rmfilenames = ["Makefile.in", "Makefile", '*.o', '*.deps/*']
for rmfile in rmfilenames:
    for filename in glob.glob(rmfile):
        os.remove(filename)

shutil.rmtree('dune')
shutil.rmtree('src')

# set CMake File for each directory
os.chdir(module_path)
def __generate_new_content(lines, pattern, replace_str):
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

def __generate_subdirectory_content(dirs):
    content = ""
    for name in dirs:
        content += "add_subdirectory(" + name + ")" + "\n"
    return content

def __generate_install_content(header_files, destination):
    if len(header_files) == 0:
        return ""
    content = "install(FILES" + "\n"
    for header_file in header_files:
        content += "    " + header_file + "\n"
    content += "DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}" + destination + ")" + "\n"
    return content

def __generate_content(cmake_list_txt_file, dirs, header_files, destination):
    subdirectory_content = __generate_subdirectory_content(dirs)
    install_content = __generate_install_content(header_files, destination)

    if os.path.exists(cmake_list_txt_file):
        with open(cmake_list_txt_file, "r", encoding="utf-8") as f:
            content = "".join(f.readlines()).strip()
            f.close()
        content = __generate_new_content(content, "add_subdirectory(", subdirectory_content)
        return __generate_new_content(content, "install(FILE", install_content)
    else:
        if subdirectory_content == "":
            return install_content
        else:
            return subdirectory_content + "\n" + install_content

def __generate_cmake_lists_txt(cmake_list_txt_file, dirs, header_files, destination):
    content = __generate_content(cmake_list_txt_file, dirs, header_files, destination)

    with open(cmake_list_txt_file, "w", encoding="utf-8") as f:
        if content != "":
            f.write(content)
        f.close()

def __check_dir(root_dir):
    if not os.path.exists(root_dir):
        print("root path" + str(root_dir) + "is not exist!")
        return False
    if not os.path.isdir(root_dir):
        print("root path" + str(root_dir) + "is not dir!")
        return False
    return True

def __check_str(root_dir):
    if root_dir is None or root_dir.strip() == "":
        print("root path is None!")
        return None
    root_dir = root_dir.strip()
    if root_dir.endswith(os.sep):
        root_dir = root_dir[:-1]
    return root_dir

def generate_cmake_lists_txt_file(root_dir):
    root_dir = __check_str(root_dir)
    if root_dir is None:
        return
    if not __check_dir(root_dir):
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
        __generate_cmake_lists_txt(cmake_list_txt_file, dirs, header_files, destination)

generate_cmake_lists_txt_file(module_path)

# move patches folder into module if existing
if (os.path.isdir("patches") ):
    subprocess.call(["mv","patches",module_name],shell=True)

# output guidence for users
print("\n"+"*"*80+"\n"
      "The extracted module is contained in the subfolder \""+module_name+"\".\n"
      "You can build it using \"dunecontrol ... --only="+ module_name +" all\".\n"
      +"*"*80+"\n"
      "BEFORE building, you can add the module to dumux-pub by something like:\n"
      "(Rename module name if it does not match the AuthorLastNameYearx scheme\n"
      "and commit it to the Git repository dumux-pub)\n"
      "git clone https://git.iws.uni-stuttgart.de/dumux-pub/AuthorLastNameYearx.git \n"
      "sed -i '/Module:/c\Module: AuthorLastNameYearx' "+ module_name+"/dune.module\n"
      "mv "+ module_name +"/* dumux-pub/AuthorLastNameYearx/.\n"
      "cd AuthorLastNameYearx\n"
      "git commit -a\n"
      "git push")

sys.exit(0)