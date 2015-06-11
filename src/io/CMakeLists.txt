#
#  Source files and headers in this directory
#
file(GLOB IO_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB IO_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${IO_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(io      SHARED ${IO_SRC})
ADD_LIBRARY(io_stat STATIC ${IO_SRC})


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(io      ${ext_libs})
TARGET_LINK_LIBRARIES(io_stat ${ext_libs})

                                                 
