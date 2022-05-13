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


target_link_libraries(libint2_wrappers PRIVATE Libint2::int2)
target_link_libraries(libint2_wrappers_stat PRIVATE Libint2::int2)

#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(libint2_wrappers      util_stat linalg_stat meigen_stat specialfunctions_stat ${ext_libs} ${Libint2_LIBRARY})
TARGET_LINK_LIBRARIES(libint2_wrappers_stat util_stat linalg_stat meigen_stat specialfunctions_stat ${ext_libs} ${Libint2_LIBRARY})


