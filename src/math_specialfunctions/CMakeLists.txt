#
#  Source files and headers in this directory
#
file(GLOB SPFUNC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB SPFUNC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${SPFUNC_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(specialfunctions      SHARED ${SPFUNC_SRC})
ADD_LIBRARY(specialfunctions_stat STATIC ${SPFUNC_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in math_specialfunctions")
   target_precompile_headers(specialfunctions REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(specialfunctions_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(specialfunctions      linalg meigen random ${ext_libs})
TARGET_LINK_LIBRARIES(specialfunctions_stat linalg_stat meigen_stat random_stat ${ext_libs})

                                                 
