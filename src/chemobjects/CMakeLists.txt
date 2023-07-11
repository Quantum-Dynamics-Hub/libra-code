#
# Set the libraries
# 
#SET( ext_libs ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}/libpython2.6.dll )


#
#  Building sub-directories
#
MESSAGE("Going into subdirectory universe...")
ADD_SUBDIRECTORY("universe")

MESSAGE("Going into subdirectory mol...")
ADD_SUBDIRECTORY("mol")

MESSAGE("Going into subdirectory chemsys...")
ADD_SUBDIRECTORY("chemsys")




#
#  Source files and headers in this directory
#
file(GLOB CHEMOBJECTS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB CHEMOBJECTS_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${CHEMOBJECTS_HEADERS}) 


#
#  Build both static and dynamic libraries
#
ADD_LIBRARY(chemobjects      SHARED ${CHEMOBJECTS_SRC})
ADD_LIBRARY(chemobjects_stat STATIC ${CHEMOBJECTS_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: chemobjects")
   target_precompile_headers(chemobjects REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(chemobjects_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(chemobjects      chemsys mol universe linalg)   
TARGET_LINK_LIBRARIES(chemobjects_stat chemsys_stat mol_stat universe_stat linalg_stat)  



