#
#  Source files and headers in this directory
#
file(GLOB MMATH_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB MMATH_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${MMATH_HEADERS}) 


#
#  External libraries: Boost and Python
#
SET( ext_libs 
     ${Boost_LIBRARIES}
     ${PYTHON_LIBRARIES}/libpython2.6.dll
   )


#
#  Build the library
#
ADD_LIBRARY(mmath SHARED ${MMATH_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(mmath ${ext_libs})

                                                 
