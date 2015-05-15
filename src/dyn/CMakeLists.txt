#
# Set the libraries
# 
#SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )


#
#  Building sub-directories
#
MESSAGE("Going into subdirectory nuclear...")
ADD_SUBDIRECTORY("nuclear")

MESSAGE("Going into subdirectory rigidbody...")
ADD_SUBDIRECTORY("rigidbody")


MESSAGE("Going into subdirectory electronic...")
ADD_SUBDIRECTORY("electronic")


#
#  Source files and headers in this directory
#
file(GLOB DYN_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB DYN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${DYN_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(dyn      SHARED ${DYN_SRC})
ADD_LIBRARY(dyn_stat STATIC ${DYN_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(dyn      electronic_stat rigidbody_stat nuclear_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(dyn_stat electronic_stat rigidbody_stat nuclear_stat mmath_stat ${ext_libs})



