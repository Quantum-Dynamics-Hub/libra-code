# Set the compiler as described here: https://cmake.org/Wiki/CMake_FAQ#How_do_I_use_a_different_compiler.3F
#
#  Verify that Cmake is up to date
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.24)  # used to be 3.12
set (CMAKE_CXX_STANDARD 17)
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
# Define USING_PCH globally for all targets
add_compile_definitions(USING_PCH)
add_compile_definitions(pch_EXPORTS)
#
option(USING_PCH "Build using pre-compiled header support" ON)
option(WITH_OPENMP "Use OpenMP for parallelization" ON)


#
#  OpenMP
#
if(WITH_OPENMP)
  find_package(OpenMP REQUIRED)
endif()


#
#  Python library
#
MESSAGE(STATUS "Looking for Python libraries...")
#FIND_PACKAGE(Python3 3.6 REQUIRED COMPONENTS Interpreter Development)
FIND_PACKAGE(Python3 3.6 REQUIRED COMPONENTS Development)

#IF(Python3_FOUND)
#  MESSAGE("Success!")
#  INCLUDE_DIRECTORIES("${PYTHON_INCLUDE_DIRS}")
IF(NOT Python3_FOUND)
  MESSAGE(STATUS "Unable to find correct Python version. If your Python include/libraries are installed in a non-standard location")
  MESSAGE(STATUS "Try using -DPython3_ROOT_DIR=\${path_to_the_python_root_dir}")
  MESSAGE(FATAL_ERROR)
ENDIF()

MESSAGE(STATUS "Found Python libraries directory: ${Python3_LIBRARIES}")
MESSAGE(STATUS "Found Python include directory:   ${Python3_INCLUDE_DIRS}")
SET(LIBRA_Python3_VERSION "${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}")




#
#  Boost library
#

MESSAGE(STATUS "Looking for Boost libraries...")
SET(Boost_NO_BOOST_CMAKE TRUE) # to fix the problem of finding python37 libs
SET(Boost_USE_STATIC_LIBS=OFF)
SET(Boost_USE_MULTITHREADED=OFF)
SET(Boost_USE_STATIC_RUNTIME=OFF)

FIND_PACKAGE(Boost 1.74.0 REQUIRED COMPONENTS regex python${LIBRA_Python3_VERSION})
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}")
ENDIF()

MESSAGE(STATUS "Found Boost libraries directory: ${Boost_LIBRARIES}")
MESSAGE(STATUS "Found Boost include directory:   ${Boost_INCLUDE_DIRS}")


#
#   Eigen 3 library
#

MESSAGE(STATUS "Looking for Eigen3 library...")
FIND_PACKAGE(Eigen3)
IF(EIGEN3_FOUND)
  INCLUDE_DIRECTORIES("${EIGEN3_INCLUDE_DIR}")
  SET(EIGEN3_INCLUDE_DIRS "${EIGEN3_INCLUDE_DIR}")
  MESSAGE(STATUS "Found Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}")
ELSEIF(NOT EIGEN3_FOUND)
  MESSAGE(FATAL_ERROR "Unable to find Eigen3")
ENDIF()


#
#  libint2 library
#

MESSAGE(STATUS "Looking for libint2 libraries...")
FIND_PACKAGE(Libint2 CONFIG 2.7.1 REQUIRED COMPONENTS shared gss impure_sh onebody_d0_l6 g12_d0_l4 g12_d1_l4 eri_c4_d0_l5 eri_c4_d1_l4) #gss e5 g5)
if(Libint2_FOUND)
  MESSAGE(STATUS "Found libint2 include directory: ${LIBINT2_INCLUDE_DIRS}")
  INCLUDE_DIRECTORIES(Libint2::int2 Libint2::cxx Libint2::int2-cxx Libint2::impure_sh Libint2::gss Libint2::onebody_d0_l6 Libint::g12_d0_l4 Libint2::g12_d1_l4 Libint2::shared Libint2::cxx_ho Libint2::c  Libint2::eri_c4_d0_l5  Libint2::eri_c4_d1_l4)
endif()


#
#  LibTorch
#

# Check your PyTorch executable:
# python -c "import torch; print(torch._C._GLIBCXX_USE_CXX11_ABI)"
# if it yields True - You have a C++11 ABI build, choose that option
# if it yields False - Your PyTorch uses the older pre-C++11 ABI, choose that option

# CXX11
SET(Torch_Root /home/alexvakimov/SOFTWARE/libpytorch/libtorch)
#ADD_DEFINITIONS(-D_GLIBCXX_USE_CXX11_ABI=1)
#export CXXFLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

# pre-CXX11 - doesn't work
#SET(Torch_Root /home/alexvakimov/SOFTWARE/libpytorch/libtorch-pre-cxx11-abi)
#ADD_DEFINITIONS(-D_GLIBCXX_USE_CXX11_ABI=0)

MESSAGE(STATUS "Looking for LibTorch libraries...")
# Set the LibTorch directory path
SET(Torch_DIR ${Torch_Root}/share/cmake/Torch)
#SET(CMAKE_PREFIX_PATH /home/alexvakimov/SOFTWARE/libpytorch/libtorch)  # Here - need to define the prefix 
#MESSAGE(${CMAKE_PREFIX_PATH})
FIND_PACKAGE(Torch REQUIRED)
SET(TORCH_LIBRARIES "${TORCH_LIBRARIES};${Torch_Root}/lib/libtorch_python.so" )
MESSAGE("Torch Libraries are: ${TORCH_LIBRARIES}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")



#
# GNU compiler definitions
#
IF(CMAKE_COMPILER_IS_GNUCXX)
# For Linux
  ADD_DEFINITIONS("-Wall -Wl,-z,defs")
  add_compile_options(-Wno-unused-variable)
  add_compile_options(-Wno-sign-compare)
# For Cygwin
#  ADD_DEFINITIONS("-Wall -DCYGWIN")
ELSE()
  MESSAGE(FATAL_ERROR "CMakeLists.txt has not been tested/written for your compiler.")
ENDIF()



#
#  Cmake configuration
#
MESSAGE(STATUS "Setting up the CMake configuration...")
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "RELEASE")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
ENDIF()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")



#
# Set the libraries
#
if(NOT WITH_OPENMP)
  SET(ext_libs Boost::regex Boost::python${LIBRA_Python3_VERSION} Python3::Python ${TORCH_LIBRARIES} torch torch_cpu )
else()
  SET(ext_libs Boost::regex Boost::python${LIBRA_Python3_VERSION} Python3::Python ${TORCH_LIBRARIES} torch torch_cpu OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()



#
# Now building the project
#
MESSAGE("Going into subdirectory src...")
ADD_SUBDIRECTORY("src")


