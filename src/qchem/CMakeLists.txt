#
# Set the libraries
# 
#SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )


#
#  Building sub-directories
#
MESSAGE("Going into subdirectory molint...")
ADD_SUBDIRECTORY("molint")

MESSAGE("Going into subdirectory qobjects...")
ADD_SUBDIRECTORY("qobjects")



#
#  Source files and headers in this directory
#
file(GLOB QCHEM_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB QCHEM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${QCHEM_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(qchem      SHARED ${QCHEM_SRC})
ADD_LIBRARY(qchem_stat STATIC ${QCHEM_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(qchem      molint_stat mmath_stat ${ext_libs})
TARGET_LINK_LIBRARIES(qchem_stat molint_stat mmath_stat ${ext_libs})



