# Import the os module, for the os.walk function
import os
import re

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
