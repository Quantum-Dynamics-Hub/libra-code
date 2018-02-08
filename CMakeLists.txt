# Set the compiler as described here: https://cmake.org/Wiki/CMake_FAQ#How_do_I_use_a_different_compiler.3F
#
#  Verify that Cmake is up to date
#
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_LEGACY_CYGWIN_WIN32 0)


#
#  Project name
#


PROJECT("libra")

#
#  Cmake configuration
#
IF(NOT CMAKE_BUILD_TYPE)
  #SET(CMAKE_BUILD_TYPE "DEBUG")
  SET(CMAKE_BUILD_TYPE "RELEASE")
  #SET(CMAKE_BUILD_TYPE "RELWITHDEBINFO")
  #SET(CMAKE_BUILD_TYPE "MINSIZEREL")
  SET(CMAKE_CXX_FLAGS "-fPIC -O2 -fpermissive -w -g")
ENDIF()

#
#  User-defined Find modules
#
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)


#
#  Boost library
#

MESSAGE("Looking for Boost libraries...")
FIND_PACKAGE(Boost 1.60.0 REQUIRED)
IF(Boost_FOUND)
  MESSAGE("Success!")
  INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}")
  SET(Boost_USE_STATIC_LIBS OFF)
  SET(Boost_USE_MULTITHREADED OFF)
  SET(Boost_USE_STATIC_RUNTIME OFF)
  FIND_PACKAGE(Boost 1.60.0 COMPONENTS python regex) # it seems that it is important to find sub-packages too!!!
ELSEIF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find correct Boost version. Did you set BOOST_ROOT?")
ENDIF()

MESSAGE("Found Boost libraries directory: ")
MESSAGE(${Boost_LIBRARIES})
MESSAGE("Found Boost include directory: ")
MESSAGE(${Boost_INCLUDE_DIRS})


#
#  Python library
#
MESSAGE("Looking for Python libraries...")
MESSAGE("Advice: If this search fails, try setting cache variables: -DPYTHON_LIBRARY=\${path_to_the_python_library} -DPYTHON_INCLUDE_DIR=\${path_to_where_Python.h_is_found}")
MESSAGE("Hint to self: -DPYTHON_INCLUDE_DIR=/home/Alexey_2/Soft/Python-2.6.7/bin/include/python2.6 -DPYTHON_LIBRARY=/home/Alexey_2/Soft/Python-2.6.7/bin/bin")

# For Linux and Cygwin
FIND_PACKAGE(PythonLibs 2.7.0 REQUIRED)

IF(PYTHONLIBS_FOUND)
  MESSAGE("Success!")
  INCLUDE_DIRECTORIES("${PYTHON_INCLUDE_DIRS}")
ELSEIF(NOT PYTHONLIBS_FOUND)
  MESSAGE("Unable to find correct Python version. If your Python include/libraries are installed in a non-standard location")
  MESSAGE("Try using -DPYTHON_LIBRARY=\${path_to_the_python_library} -PYTHON_INCLUDE_DIR=\${path_to_where_Python.h_is_found}")
  MESSAGE(FATAL_ERROR)
ENDIF()

MESSAGE("Found  Python libraries directory: ")
MESSAGE(${PYTHON_LIBRARIES})
MESSAGE("Found Python include directory: ")
MESSAGE(${PYTHON_INCLUDE_DIRS})


#
# GNU compiler definitions
#
IF(CMAKE_COMPILER_IS_GNUCXX)
# For Linux
#  ADD_DEFINITIONS("-Wall")
# For Cygwin
  ADD_DEFINITIONS("-Wall -DCYGWIN")
ELSE()
  MESSAGE(FATAL_ERROR "CMakeLists.txt has not been tested/written for your compiler.")
ENDIF()



#
# Set the libraries
# 
SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ) 



#
# Now building the project
#
MESSAGE("Going into subdirectory src...")
ADD_SUBDIRECTORY("src")


