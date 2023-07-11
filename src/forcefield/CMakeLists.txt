#
#  Source files and headers in this directory
#
file(GLOB FORCE_FIELD_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB FORCE_FIELD_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${FORCE_FIELD_HEADERS}) 


#SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(forcefield SHARED ${FORCE_FIELD_SRC})
ADD_LIBRARY(forcefield_stat STATIC ${FORCE_FIELD_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: forcefield")
   target_precompile_headers(forcefield REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(forcefield_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#

TARGET_LINK_LIBRARIES( forcefield      io linalg specialfunctions)
TARGET_LINK_LIBRARIES( forcefield_stat io_stat linalg_stat specialfunctions_stat)


