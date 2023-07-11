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



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: calculators")
   target_precompile_headers(calculators REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(calculators_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(calculators  linalg meigen specialfunctions common_types timer)

TARGET_LINK_LIBRARIES(calculators_stat linalg_stat meigen_stat specialfunctions_stat  common_types_stat timer_stat)

                                                 
