#
#  Source files and headers in this directory
#
file(GLOB UNIVERSE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB UNIVERSE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${UNIVERSE_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(universe SHARED ${UNIVERSE_SRC})
ADD_LIBRARY(universe_stat STATIC ${UNIVERSE_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: universe")
   target_precompile_headers(universe REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(universe_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(universe      io)
TARGET_LINK_LIBRARIES(universe_stat io_stat)


