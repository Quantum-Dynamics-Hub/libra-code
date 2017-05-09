#
#  Source files and headers in this directory
#
file(GLOB CALCULATORS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB CALCULATORS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CALCULATORS_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(calculators      SHARED ${CALCULATORS_SRC})
ADD_LIBRARY(calculators_stat STATIC ${CALCULATORS_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(calculators      linalg_stat meigen_stat specialfunctions_stat
                      common_types_stat timer_stat ${ext_libs})

TARGET_LINK_LIBRARIES(calculators_stat linalg_stat meigen_stat specialfunctions_stat
                      common_types_stat timer_stat ${ext_libs})

                                                 
