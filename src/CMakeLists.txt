#
# Set the libraries
# 
SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )

#SET( ext_libs1 ${CMAKE_CURRENT_SOURCE_DIR}/cygmmath.dll ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )
#
#  mmath
#
#file(GLOB MMATH_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/mmath/*.h)
#file(GLOB MMATH_SRC ${CMAKE_CURRENT_SOURCE_DIR}/mmath/*.cpp ${MMATH_HEADERS}) 
#ADD_LIBRARY(mmath STATIC ${MMATH_SRC})
#
#
#  ann
#
#file(GLOB ANN_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/ann/*.h)
#file(GLOB ANN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/ann/*.cpp ${ANN_HEADERS}) 
#ADD_LIBRARY(ann STATIC ${ANN_SRC})
#


#
#  Set the project-wide dependencies:
#
#TARGET_LINK_LIBRARIES(mmath ${ext_libs})
#TARGET_LINK_LIBRARIES(ann ${ext_libs1} )




#
#  Building sub-directories
#
MESSAGE("Going into subdirectory mmath...")
ADD_SUBDIRECTORY("mmath")

MESSAGE("Going into subdirectory ann...")
ADD_SUBDIRECTORY("ann")

