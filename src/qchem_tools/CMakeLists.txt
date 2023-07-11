#
# Set the libraries
# 
#SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )


#
#  Source files and headers in this directory
#
file(GLOB QCHEM_TOOLS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB QCHEM_TOOLS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${QCHEM_TOOLS_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(qchem_tools      SHARED ${QCHEM_TOOLS_SRC})
ADD_LIBRARY(qchem_tools_stat STATIC ${QCHEM_TOOLS_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: qchem_tools")
   target_precompile_headers(qchem_tools REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(qchem_tools_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()





#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(qchem_tools      chemobjects
                      atomistic qobjects
                      basis_setups control_parameters
                      model_parameters
                      linalg ${ext_libs} )

TARGET_LINK_LIBRARIES(qchem_tools_stat      chemobjects_stat 
                      atomistic_stat qobjects_stat
                      basis_setups_stat control_parameters_stat
                      model_parameters_stat
                      linalg_stat ${ext_libs} )



