#
#  Source files and headers in this directory
#
file(GLOB ATOMISTIC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB ATOMISTIC_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${ATOMISTIC_HEADERS}) 


MESSAGE("Going into subdirectory Hamiltonian_Generic...")
ADD_SUBDIRECTORY("Hamiltonian_Generic")

MESSAGE("Going into subdirectory Hamiltonian_MM...")
ADD_SUBDIRECTORY("Hamiltonian_MM")


MESSAGE("Going into subdirectory Hamiltonian_QM...")
ADD_SUBDIRECTORY("Hamiltonian_QM")


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(atomistic SHARED ${ATOMISTIC_SRC})
ADD_LIBRARY(atomistic_stat STATIC ${ATOMISTIC_SRC})



if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: atomistic")
   target_precompile_headers(atomistic REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(atomistic_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()



#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES( atomistic 
                       common_types
                       calculators
                       hamiltonian_mm
                       hamiltonian_qm
                       hamiltonian_generic
                       chemobjects
                       linalg meigen)


TARGET_LINK_LIBRARIES( atomistic_stat
                       common_types_stat
                       calculators_stat
                       hamiltonian_mm_stat 
                       hamiltonian_qm_stat
                       hamiltonian_generic_stat
                       chemobjects_stat
                       linalg_stat meigen_stat)

