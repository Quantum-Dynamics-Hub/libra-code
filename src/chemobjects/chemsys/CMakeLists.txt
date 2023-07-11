#
#  Source files and headers in this directory
#
file(GLOB CHEMSYS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB CHEMSYS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CHEMSYS_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(chemsys SHARED ${CHEMSYS_SRC})
ADD_LIBRARY(chemsys_stat STATIC ${CHEMSYS_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: chemsys")
   target_precompile_headers(chemsys REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(chemsys_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(chemsys      mol linalg random graph)
TARGET_LINK_LIBRARIES(chemsys_stat mol_stat linalg_stat random_stat graph_stat)


