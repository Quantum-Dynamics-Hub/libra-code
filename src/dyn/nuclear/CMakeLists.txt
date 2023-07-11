#
#  Source files and headers in this directory
#
file(GLOB NUCLEAR_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB NUCLEAR_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${NUCLEAR_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(nuclear SHARED ${NUCLEAR_SRC})
ADD_LIBRARY(nuclear_stat STATIC ${NUCLEAR_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in nuclear")
   target_precompile_headers(nuclear REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(nuclear_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(nuclear      linalg_stat ${ext_libs})
TARGET_LINK_LIBRARIES(nuclear_stat linalg_stat ${ext_libs})


