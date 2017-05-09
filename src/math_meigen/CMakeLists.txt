#
#
#
#
#MESSAGE("Going into subdirectory Eigen/eigen-eigen-b9cd8366d4e8...")
#ADD_SUBDIRECTORY("Eigen/eigen-eigen-b9cd8366d4e8")

#MESSAGE("Going into subdirectory Eigen/eigen-eigen-b9cd8366d4e8...")
#ADD_SUBDIRECTORY("Eigen")


# Make the execution of this script conditional, depending on the Eigen library
# availablility
# In general, we need to call the find_Eigen module
# but for now, I just set the include dir manually
#
#set (EIGEN_INCLUDE /Eigen/eigen-eigen-b9cd8366d4e8/ )
set (EIGEN_INCLUDE ./ )

#
#  Source files and headers in this directory
#
file(GLOB EIG_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB EIG_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${EIG_HEADERS}) 


INCLUDE_DIRECTORIES(${EIGEN_INCLUDE})




#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(meigen SHARED ${EIG_SRC})
ADD_LIBRARY(meigen_stat STATIC ${EIG_SRC})



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(meigen      linalg_stat ${ext_libs})
TARGET_LINK_LIBRARIES(meigen_stat linalg_stat ${ext_libs})


