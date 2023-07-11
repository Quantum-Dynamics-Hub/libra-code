#
#  Source files and headers in this directory
#
file(GLOB LIBINT2_WRAPPERS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB LIBINT2_WRAPPERS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${LIBINT2_WRAPPERS_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(libint2_wrappers SHARED ${LIBINT2_WRAPPERS_SRC})
ADD_LIBRARY(libint2_wrappers_stat STATIC ${LIBINT2_WRAPPERS_SRC})



#
#  PCH (Precompiled header) support
# 
option(USING_PCH "Build using pre-compiled header support" ON)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support")
   target_precompile_headers(libint2_wrappers REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(libint2_wrappers_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
   set(CMAKE_PCH_INSTANTIATE_TEMPLATES ON)
   add_compile_definitions(USING_PCH)
else()
   message(STATUS "No pre-compiled header support")
endif()




target_link_libraries(libint2_wrappers PRIVATE Libint2::int2)
target_link_libraries(libint2_wrappers_stat PRIVATE Libint2::int2)

#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(libint2_wrappers      util_stat linalg_stat meigen_stat specialfunctions_stat ${ext_libs} ${Libint2_LIBRARY})
TARGET_LINK_LIBRARIES(libint2_wrappers_stat util_stat linalg_stat meigen_stat specialfunctions_stat ${ext_libs} ${Libint2_LIBRARY})


