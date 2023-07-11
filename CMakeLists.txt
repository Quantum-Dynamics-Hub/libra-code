# Set the compiler as described here: https://cmake.org/Wiki/CMake_FAQ#How_do_I_use_a_different_compiler.3F
#
#  Verify that Cmake is up to date
#
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
#SET(CMAKE_LEGACY_CYGWIN_WIN32 0)


#
#  Project name
#


project("libra")


#
#  User-defined Find modules
#
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)


#
# Try using the PCH - it will only work if you have the right version of cmake
#
option(USING_PCH "Build using pre-compiled header support" ON)


#
#  Python library
#
MESSAGE("Looking for Python libraries...")
MESSAGE("Advice: If this search fails, try setting cache variables: -DPYTHON_LIBRARY=\${path_to_the_python_library} -DPYTHON_INCLUDE_DIR=\${path_to_where_Python.h_is_found}")
MESSAGE("Hint to self: -DPYTHON_INCLUDE_DIR=/home/Alexey_2/Soft/Python-2.6.7/bin/include/python2.6 -DPYTHON_LIBRARY=/home/Alexey_2/Soft/Python-2.6.7/bin/bin")

# For Linux and Cygwin
# Adapt it according to your system
#FIND_PACKAGE(PythonLibs 3.6 REQUIRED)
FIND_PACKAGE(PythonLibs 3.7 REQUIRED)
#FIND_PACKAGE(PythonLibs 3.7 REQUIRED)

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
#  Boost library
#

MESSAGE("Looking for Boost libraries...")
SET(Boost_NO_BOOST_CMAKE TRUE) # to fix the problem of finding python37 libs
SET(Boost_USE_STATIC_LIBS=OFF)
SET(Boost_USE_MULTITHREADED=OFF)
SET(Boost_USE_STATIC_RUNTIME=OFF)

# Adapt it according to your system
#FIND_PACKAGE(Boost 1.73.0 REQUIRED)
FIND_PACKAGE(Boost 1.80.0 REQUIRED)
#FIND_PACKAGE(Boost REQUIRED)
IF(Boost_FOUND)
  MESSAGE("Success!")
  INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}") 
  # Adapt it based on your system
  #FIND_PACKAGE(Boost COMPONENTS python36 regex) # it seems that it is important to find sub-packages too!!!
  FIND_PACKAGE(Boost COMPONENTS python37 regex)
  #FIND_PACKAGE(Boost COMPONENTS python37 regex)
ELSEIF(NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find correct Boost version. Did you set BOOST_ROOT?")
ENDIF()

MESSAGE("Found Boost libraries directory: ")
MESSAGE("${Boost_LIBRARIES}")
MESSAGE("Found Boost include directory: ")
MESSAGE("${Boost_INCLUDE_DIRS}")



FIND_PACKAGE(Eigen3)
IF(EIGEN3_FOUND)
  INCLUDE_DIRECTORIES("${EIGEN3_INCLUDE_DIR}")   
  MESSAGE("Found Eigen3 include directory: ")
  MESSAGE("${EIGEN3_INCLUDE_DIR}")
ELSEIF(NOT EIGEN3_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find Eigen3")
ENDIF()


SET(EIGEN3_INCLUDE_DIRS "${EIGEN3_INCLUDE_DIR}")

# Remove CMakeCache.txt from the previous build, if you get an error finding libint2 due to problems 
# in finding Eigen3

# Libint2 2.7.0 CONFIG REQUIRED COMPONENTS shared gss e5 g5
#FIND_PACKAGE(Libint2 MODULE 2.7.0)
#find_package(Libint2)
FIND_PACKAGE(Libint2 CONFIG 2.7.1 REQUIRED COMPONENTS shared gss impure_sh onebody_d0_l6 g12_d0_l4 g12_d1_l4 eri_c4_d0_l5 eri_c4_d1_l4) #gss e5 g5)
IF(Libint2_FOUND)
  #INCLUDE_DIRECTORIES("${Libint2_INCLUDE_DIRS}")   
  INCLUDE_DIRECTORIES(Libint2::int2 Libint2::cxx Libint2::int2-cxx Libint2::impure_sh Libint2::gss Libint2::onebody_d0_l6 Libint::g12_d0_l4 Libint2::g12_d1_l4 Libint2::shared Libint2::cxx_ho Libint2::c  Libint2::eri_c4_d0_l5  Libint2::eri_c4_d1_l4)   
  MESSAGE("${LIBINT2_INCLUDE_DIR}")
  MESSAGE("Found Libint2 include directory: ")
  MESSAGE("${Libint2_INCLUDE_DIRS}")
  MESSAGE("Found Libint2_LIBRARIES: ")
  MESSAGE("${Libint2_LIBRARIES}")
  MESSAGE("Found Libint2_LIBRARY: ")
  MESSAGE("${Libint2_LIBRARY}")
ELSEIF(NOT Libint2_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find correct libint2")
ENDIF()


#
# GNU compiler definitions
#
IF(CMAKE_COMPILER_IS_GNUCXX)
# For Linux
  ADD_DEFINITIONS("-Wall -Wl,-z,defs")
# For Cygwin
#  ADD_DEFINITIONS("-Wall -DCYGWIN")
ELSE()
  MESSAGE(FATAL_ERROR "CMakeLists.txt has not been tested/written for your compiler.")
ENDIF()



#
#  Cmake configuration
#
MESSAGE("Setting up the CMake configuration...")
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "RELEASE")
  SET(CMAKE_CXX_FLAGS "-fPIC -O2 -fpermissive -w -I ${Boost_INCLUDE_DIRS} -fopenmp")
ENDIF()



#
# Set the libraries
# 
SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}) 



#
# Now building the project
#
MESSAGE("Going into subdirectory src...")
ADD_SUBDIRECTORY("src")


