import os
import sys
import pprint


def modify2(filename):
    file = open(filename, 'r')
    file_data = file.readlines()
    file.close()
    tmp = file_data[2]
    file_data[2] = file_data[4]
    file_data[4] = tmp
    re_file = open(filename, 'w+')
    for line in file_data:
        re_file.write(line)
    re_file.close()


def edit_file(filename):
    file = open(filename, 'r')
    file_data = file.readlines()
    file.close()
    for i in range(26):
        file_data.pop(0)
    new_copyright = ["// @HEADER\n",
                     "// ****************************************************************************\n",
                     "//                Hydrofem: Copyright (2016) S. Mabuza\n",
                     "//\n",
                     "// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)\n",
                     "// ****************************************************************************\n",
                     "// @HEADER\n"]
    re_file = open(filename, 'w+')
    for line in new_copyright:
        re_file.write(line)
    for line in file_data:
        re_file.write(line)
    re_file.close()

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # print(sys.argv)
    edit_file(sys.argv[1])
    # modify2(sys.argv[1])

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
