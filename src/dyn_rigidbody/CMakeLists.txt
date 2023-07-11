#
#  Source files and headers in this directory
#
file(GLOB RB_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB RB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${RB_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(rigidbody SHARED ${RB_SRC})
ADD_LIBRARY(rigidbody_stat STATIC ${RB_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in rigidbody")
   target_precompile_headers(rigidbody REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(rigidbody_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(rigidbody      linalg meigen specialfunctions io)
TARGET_LINK_LIBRARIES(rigidbody_stat linalg_stat meigen_stat specialfunctions_stat io_stat)


