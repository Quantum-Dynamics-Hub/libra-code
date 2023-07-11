#
# Set the libraries
# 
#SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )


#
#  Building sub-directories
#
MESSAGE("Going into subdirectory state...")
ADD_SUBDIRECTORY("state")


#
#  Source files and headers in this directory
#
file(GLOB SRIPTS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB SRIPTS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${SRIPTS_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(scripts      SHARED ${SRIPTS_SRC})
ADD_LIBRARY(scripts_stat STATIC ${SRIPTS_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(scripts      state ${ext_libs})
TARGET_LINK_LIBRARIES(scripts_stat state_stat ${ext_libs})



