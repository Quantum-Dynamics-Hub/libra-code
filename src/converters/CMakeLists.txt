#
#  Source files and headers in this directory
#
file(GLOB CONVERTERS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB CONVERTERS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CONVERTERS_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(converters      SHARED ${CONVERTERS_SRC})
ADD_LIBRARY(converters_stat STATIC ${CONVERTERS_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in converters")
   target_precompile_headers(converters REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(converters_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(converters      dyn chemobjects linalg)
TARGET_LINK_LIBRARIES(converters_stat dyn_stat chemobjects_stat linalg_stat)



