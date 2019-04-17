#*********************************************************************************
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file path_libra_lib.py
# This module defines the function which makes paths to the libra libraries.
# (If you can import {cyg,lib}libra_core, you don't need to use this.)

import os

def path_libra_lib(librapath):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] librapath  the path to the source directory containing libra libraries.
    # This returns the paths to libra modules using environmetal variables.
    #
    # Used in: run.py

    # detect all directories
    libdirs = []
    for d in os.listdir(librapath):
        if os.path.isdir(os.path.join(librapath,d)):
            libdirs.append(d)

    for d in libdirs:
        if not d == "CMakeFiles": # need not path "CMakefiles" directory.
            os.environ["libra_%s_path"%(d)] = librapath + "/" + d
            os.system("echo libra_%s_path = $libra_%s_path" %(d,d))
            #print d

