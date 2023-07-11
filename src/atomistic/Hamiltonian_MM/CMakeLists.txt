#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_MM_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_MM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_MM_HEADERS}) 


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian_mm SHARED ${HAMILTONIAN_MM_SRC})
ADD_LIBRARY(hamiltonian_mm_stat STATIC ${HAMILTONIAN_MM_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in hamiltonian_mm")
   target_precompile_headers(hamiltonian_mm REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(hamiltonian_mm_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( hamiltonian_mm  
                       chemobjects
                       forcefield
                       pot
                       cell
                       linalg specialfunctions meigen)

TARGET_LINK_LIBRARIES( hamiltonian_mm_stat 
                       chemobjects_stat
                       forcefield_stat
                       pot_stat
                       cell_stat
                       linalg_stat specialfunctions_stat meigen_stat)


