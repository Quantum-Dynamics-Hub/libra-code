#
#  Source files and headers in this directory
#
file(GLOB IVR_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB IVR_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${IVR_HEADERS}) 

#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(ivr SHARED ${IVR_SRC})
ADD_LIBRARY(ivr_stat STATIC ${IVR_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(ivr      meigen_stat random_stat linalg_stat specialfunctions_stat ${ext_libs})
TARGET_LINK_LIBRARIES(ivr_stat meigen_stat random_stat linalg_stat specialfunctions_stat ${ext_libs})


