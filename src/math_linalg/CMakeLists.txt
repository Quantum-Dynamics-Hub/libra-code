#
#  Source files and headers in this directory
#
file(GLOB LINALG_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB LINALG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${LINALG_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(linalg      SHARED ${LINALG_SRC})
ADD_LIBRARY(linalg_stat STATIC ${LINALG_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in math_linalg")
   target_precompile_headers(linalg REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(linalg_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()

#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(linalg      io ${ext_libs})
TARGET_LINK_LIBRARIES(linalg_stat io_stat ${ext_libs})

                                                 
