
#
#  Source files and headers in this directory
#
file(GLOB HAMILTONIAN_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB HAMILTONIAN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ${HAMILTONIAN_HEADERS}) 



MESSAGE("Going into subdirectory nHamiltonian_Generic...")
ADD_SUBDIRECTORY("nHamiltonian_Generic")

#MESSAGE("Going into subdirectory Hamiltonian_Model...")
#ADD_SUBDIRECTORY("Hamiltonian_Model")


#MESSAGE("Going into subdirectory Hamiltonian_Extern...")
#ADD_SUBDIRECTORY("Hamiltonian_Extern")

#MESSAGE("Going into subdirectory Hamiltonian_Atomistic...")
#ADD_SUBDIRECTORY("Hamiltonian_Atomistic")


#
#  Create both static and dynamic libraries
#
ADD_LIBRARY(hamiltonian SHARED ${HAMILTONIAN_SRC})
ADD_LIBRARY(hamiltonian_stat STATIC ${HAMILTONIAN_SRC})


if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16 AND USING_PCH)  # Support for PCHs in CMake was added in 3.16
   message(STATUS "Compiling using pre-compiled header support: in hamiltonian")
   target_precompile_headers(hamiltonian REUSE_FROM pch) # With PUBLIC they will be used by targets using this target
   target_precompile_headers(hamiltonian_stat REUSE_FROM pch_stat) # With PUBLIC they will be used by targets using this target
endif()


#
#  Link to external libraries
#
TARGET_LINK_LIBRARIES(hamiltonian  nhamiltonian_generic_stat  linalg_stat ${ext_libs}     )

TARGET_LINK_LIBRARIES(hamiltonian_stat   nhamiltonian_generic_stat  linalg_stat ${ext_libs}  )
