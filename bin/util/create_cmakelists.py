# Create files CMakeLists.txt for the current folder and all subfolders.
# Should be run from the folder `dumux` that contains the header files, namely,
# one level below the top dumux folder:
# python ../bin/create_cmakelists.py.

# Import the os module, for the os.walk function
import os

# Set the directory you want to start from
rootDir = '.'
for folderName, subFolders, files in os.walk(rootDir):
    subFolders = sorted(subFolders)
    files = sorted(files)

    cmakelists = open(folderName + "/CMakeLists.txt", "w")

    for subFolder in subFolders:
        cmakelists.write("add_subdirectory(\"%s\")\n" % subFolder)

    headersExist = False
    for fileName in files:
        if fileName != "CMakeLists.txt":
            headersExist = True
            break

    if headersExist:
        if subFolders:
            cmakelists.write("\n")

        cmakelists.write("install(FILES\n")

        for fileName in files:
            if fileName != "CMakeLists.txt":
                cmakelists.write("%s\n" % fileName)

        cmakelists.write("DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/%s)\n" % folderName[2:])

    cmakelists.close()
