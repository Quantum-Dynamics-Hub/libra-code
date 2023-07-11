#
#  Source files and headers in this directory
#
file(GLOB QTAG_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB QTAG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${QTAG_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(qtag SHARED ${QTAG_SRC})
ADD_LIBRARY(qtag_stat STATIC ${QTAG_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in qtag")
   target_precompile_headers(qtag REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(qtag_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(qtag      linalg_stat ${ext_libs})
TARGET_LINK_LIBRARIES(qtag_stat linalg_stat ${ext_libs})


