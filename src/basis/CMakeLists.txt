#
#  Source files and headers in this directory
#
file(GLOB BASIS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB BASIS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${BASIS_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(basis SHARED ${BASIS_SRC})
ADD_LIBRARY(basis_stat STATIC ${BASIS_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: basis")
   target_precompile_headers(basis REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(basis_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(basis      qobjects molint linalg)
TARGET_LINK_LIBRARIES(basis_stat qobjects_stat molint_stat linalg_stat)


